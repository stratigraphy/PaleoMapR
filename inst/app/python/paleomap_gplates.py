"""
Paleogeographic map editor - thin pyGPlates bridge, called from R through reticulate.

Only the plate-rotation work lives here (R has no offline equivalent of pyGPlates).
Everything else (records, DEM correction, maps, statistics) is in R.

All functions take/return plain numpy arrays or lists so reticulate converts them
to R vectors / lists automatically.

Requirements (installed by PaleoMapR::pm_setup()):  pygplates>=1.0  plate-model-manager  numpy
"""
import numpy as np
import pygplates
from plate_model_manager import PlateModelManager

# model name -> anchor plate id. Torsvik & Cocks 2017: anchor 1 = palaeomagnetic
# (PM-HYBRID) frame; anchor 0 = TPW-corrected frame.
ANCHORS = {"scotese_and_wright2018": 0, "merdith2021": 0, "cao2024": 0,
           "torsvikcocks2017": 1, "golonka": 0, "muller2022": 0, "muller2025": 0,
           "paleomap": 0}
_models = {}
_model_dir = "models"


def _out(d):
    """numpy -> plain Python lists (converts cleanly with every reticulate version)"""
    return {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in d.items()}


def init(model_dir="models"):
    """set the folder where plate-model-manager stores the plate models"""
    global _model_dir
    _model_dir = model_dir
    return True


def _load(name):
    """rotation model + static polygons of a plate model (downloaded on first use)"""
    if name not in _models:
        base = name.replace("_TPW", "")
        anchor = 0 if name.endswith("_TPW") else ANCHORS.get(base, 0)
        m = PlateModelManager().get_model(base, data_dir=_model_dir)
        rot = pygplates.RotationModel(m.get_rotation_model(), default_anchor_plate_id=anchor)
        sp = pygplates.FeatureCollection()
        for f in m.get_layer("StaticPolygons"):
            sp.add(pygplates.FeatureCollection(f))
        _models[name] = (rot, sp)
    return _models[name]


def _rotation(rot, age_to, pid, age_from=0.0):
    return rot.get_rotation(float(age_to), int(pid), float(age_from),
                            use_identity_for_missing_plate_ids=False)


def _rotate_groups(lat, lon, pid, rot, age_to, age_from):
    """rotate points grouped by plate id (fast: one MultiPoint per plate)"""
    n = len(lat)
    out_la = np.full(n, np.nan); out_lo = np.full(n, np.nan)
    for p in np.unique(pid[pid > 0]):
        R = _rotation(rot, age_to, p, age_from)
        if R is None:
            continue
        idx = np.where(pid == p)[0]
        mp = pygplates.MultiPointOnSphere(list(zip(lat[idx], lon[idx])))
        ll = (R * mp).to_lat_lon_array()
        out_la[idx] = ll[:, 0]; out_lo[idx] = ll[:, 1]
    return out_la, out_lo


def _plate_ids(partitioner, lat, lon):
    """plate id of every point (0 = outside all polygons)"""
    pid = np.zeros(len(lat), dtype=np.int64)
    for k, (la, lo) in enumerate(zip(lat, lon)):
        r = partitioner.partition_point(pygplates.PointOnSphere(float(la), float(lo)))
        if r is not None:
            pid[k] = r.get_feature().get_reconstruction_plate_id()
    return pid


def _plate_ids_fast(partitioner, lat, lon, chunk=20000):
    """plate ids for many points: partition MultiPoints and map pieces back to the
    input points by their (rounded) coordinates; leftovers fall back to point tests"""
    n = len(lat)
    pid = np.full(n, -1, dtype=np.int64)
    key = np.round(lat, 6) * 1e7 + np.round(np.mod(lon + 180.0, 360.0), 6)
    order = np.argsort(key); skey = key[order]
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        mp = pygplates.MultiPointOnSphere(list(zip(lat[s:e], lon[s:e])))
        part, unpart = [], []
        partitioner.partition_geometry(mp, part, unpart)
        for poly, pieces in part:
            p = poly.get_feature().get_reconstruction_plate_id()
            for g in pieces:
                ll = g.to_lat_lon_array()
                k = np.round(ll[:, 0], 6) * 1e7 + np.round(np.mod(ll[:, 1] + 180.0, 360.0), 6)
                pos = np.searchsorted(skey, k)
                pos = np.clip(pos, 0, n - 1)
                ok = np.abs(skey[pos] - k) < 1e-3
                pid[order[pos[ok]]] = p
        for g in unpart:
            ll = g.to_lat_lon_array()
            k = np.round(ll[:, 0], 6) * 1e7 + np.round(np.mod(ll[:, 1] + 180.0, 360.0), 6)
            pos = np.clip(np.searchsorted(skey, k), 0, n - 1)
            ok = np.abs(skey[pos] - k) < 1e-3
            sel = order[pos[ok]]
            pid[sel[pid[sel] < 0]] = 0
    left = np.where(pid < 0)[0]
    if len(left):
        pid[left] = _plate_ids(partitioner, lat[left], lon[left])
    return pid


def reconstruct(lat, lon, ages, model="scotese_and_wright2018"):
    """present-day points -> palaeo-positions at their own ages.
    Returns dict(lat, lon, plate_id); NaN where the point is on no plate or the
    plate has no rotation at that age."""
    lat = np.atleast_1d(np.asarray(lat, float)); lon = np.atleast_1d(np.asarray(lon, float))
    ages = np.broadcast_to(np.atleast_1d(np.asarray(ages, float)), lat.shape).copy()
    rot, sp = _load(model)
    part = pygplates.PlatePartitioner(sp, rot, 0.0)
    pid = _plate_ids(part, lat, lon)
    out_la = np.full(len(lat), np.nan); out_lo = np.full(len(lat), np.nan)
    for a in np.unique(ages):
        m = ages == a
        la, lo = _rotate_groups(lat[m], lon[m], pid[m], rot, a, 0.0)
        out_la[m] = la; out_lo[m] = lo
    return _out({"lat": out_la, "lon": out_lo, "plate_id": pid})


def to_present(plat, plon, age, model="scotese_and_wright2018"):
    """palaeo-positions at `age` -> present-day positions (uses the static
    polygons reconstructed to `age`). NaN if the point is on no plate at `age`."""
    plat = np.atleast_1d(np.asarray(plat, float)); plon = np.atleast_1d(np.asarray(plon, float))
    rot, sp = _load(model)
    rec = []
    pygplates.reconstruct(sp, rot, rec, float(age), anchor_plate_id=rot.get_default_anchor_plate_id())
    part = pygplates.PlatePartitioner(rec, rot)
    pid = _plate_ids(part, plat, plon)
    la, lo = _rotate_groups(plat, plon, pid, rot, 0.0, float(age))
    return _out({"lat": la, "lon": lo, "plate_id": pid})


def transfer(plat, plon, age, source="scotese_and_wright2018", target="torsvikcocks2017",
             valid_only=True):
    """palaeo-positions at `age` in the source model -> palaeo-positions at `age`
    in the target model (via present day). Used to move PALEOMAP topography into
    other plate models. valid_only: target polygons must exist at `age`.
    Returns dict(lat, lon) with NaN for points that cannot be moved."""
    plat = np.atleast_1d(np.asarray(plat, float)); plon = np.atleast_1d(np.asarray(plon, float))
    rot0, sp0 = _load(source)
    rec = []
    pygplates.reconstruct(sp0, rot0, rec, float(age), anchor_plate_id=rot0.get_default_anchor_plate_id())
    pid0 = _plate_ids_fast(pygplates.PlatePartitioner(rec, rot0), plat, plon)
    pla, plo = _rotate_groups(plat, plon, pid0, rot0, 0.0, float(age))
    if target == source:
        ok = np.isfinite(pla)
        return _out({"lat": np.where(ok, plat, np.nan), "lon": np.where(ok, plon, np.nan)})
    rotT, spT = _load(target)
    feats = [f for f in spT if (not valid_only) or f.is_valid_at_time(float(age))]
    partT = pygplates.PlatePartitioner(feats, rotT, 0.0)
    ok = np.isfinite(pla)
    pidT = np.zeros(len(plat), dtype=np.int64)
    if ok.any():
        pidT[ok] = _plate_ids_fast(partT, pla[ok], plo[ok])
    la, lo = _rotate_groups(np.nan_to_num(pla), np.nan_to_num(plo), pidT, rotT, float(age), 0.0)
    return _out({"lat": la, "lon": lo})


# ------------------------------------------------------------------------------
# present-day lines (Natural Earth borders / coastlines) moved to a palaeo-age
# ------------------------------------------------------------------------------
_partitioned = {}
_polys = {}


def _cap(geom):
    """bounding cap of a geometry: unit centre vector and angular radius"""
    xyz = np.array([pt.to_xyz() for pt in geom.get_points()])
    c = xyz.mean(0); nc = np.linalg.norm(c)
    c = c / nc if nc > 1e-9 else xyz[0]
    r = np.arccos(np.clip(xyz @ c, -1, 1)).max()
    return c, r


def _plate_polygons(model):
    if model not in _polys:
        rot, sp = _load(model); out = []
        for f in sp:
            pid = f.get_reconstruction_plate_id(); begin = f.get_valid_time()[0]
            for g in f.get_geometries():
                if isinstance(g, pygplates.PolygonOnSphere) and pid != 0:
                    out.append((g, pid, begin))
        caps = [_cap(g) for g, _, _ in out]
        _polys[model] = (out, np.array([c for c, _ in caps]), np.array([r for _, r in caps]))
    return _polys[model]


def _nearest_plate(model, geom, age, max_km, need_rotation):
    """plate of the nearest static polygon within max_km (bounding caps prune the search)"""
    rot, _ = _load(model)
    polys, centres, radii = _plate_polygons(model)
    thr = max_km / 6371.0
    c, r = _cap(geom)
    lower = np.arccos(np.clip(centres @ c, -1, 1)) - radii - r
    best = (None, None, thr)
    for k in np.argsort(lower):
        if lower[k] > best[2]:
            break
        poly, pid, begin = polys[k]
        if need_rotation and _rotation(rot, age, pid) is None:
            continue
        d = pygplates.GeometryOnSphere.distance(geom, poly, distance_threshold_radians=best[2])
        if d is not None and d < best[2]:
            best = (pid, begin, d)
    return best[0], best[1]


def modern_lines(shapefile, model, age, snap_km=150.0, snap_km_norot=400.0):
    """Present-day lines from a shapefile, cut into the model's plates and rotated
    to `age`. Pieces outside every polygon join the nearest plate within snap_km
    (plate polygons follow the coast, so half a coastline can fall just outside);
    pieces on plates with no rotation at `age` join the nearest plate that has one
    within snap_km_norot and are flagged 'young' (approximate position).
    Returns a list of dicts: lat, lon, young (True = crust added after `age` or
    approximate position: draw dotted)."""
    key = (shapefile, model)
    rot, sp = _load(model)
    if key not in _partitioned:
        parts = pygplates.partition_into_plates(
            sp, rot, pygplates.FeatureCollection(shapefile),
            properties_to_copy=[pygplates.PartitionProperty.reconstruction_plate_id,
                                pygplates.PartitionProperty.valid_time_period])
        pieces = []
        for ft in parts:
            pid = ft.get_reconstruction_plate_id(); begin = ft.get_valid_time()[0]
            for g in ft.get_geometries():
                if pid == 0:
                    pid2, b2 = _nearest_plate(model, g, 0.0, snap_km, False)
                    if pid2 is None:
                        continue
                    pieces.append((g, pid2, b2))
                else:
                    pieces.append((g, pid, begin))
        _partitioned[key] = pieces
    out = []
    for g, pid, begin in _partitioned[key]:
        R = _rotation(rot, age, pid); young = False
        if R is None:
            pid2, _ = _nearest_plate(model, g, age, snap_km_norot, True)
            if pid2 is None:
                continue
            R = _rotation(rot, age, pid2); young = True
        young = young or (np.isfinite(begin) and begin < age)
        ll = (R * g).to_lat_lon_array()
        out.append({"lat": ll[:, 0].tolist(), "lon": ll[:, 1].tolist(), "young": bool(young)})
    return out


def plate_info(lat, lon, model="scotese_and_wright2018"):
    """present-day point -> plate id and polygon name (for site checks)"""
    rot, sp = _load(model)
    part = pygplates.PlatePartitioner(sp, rot, 0.0)
    pid, name, begin = [], [], []
    for la, lo in zip(np.atleast_1d(lat), np.atleast_1d(lon)):
        r = part.partition_point(pygplates.PointOnSphere(float(la), float(lo)))
        f = r.get_feature() if r is not None else None
        pid.append(f.get_reconstruction_plate_id() if f else -1)
        name.append(f.get_name() if f else "")
        begin.append(float(f.get_valid_time()[0]) if f else np.nan)
    return {"plate_id": [int(x) for x in pid], "name": name, "valid_from": begin}
