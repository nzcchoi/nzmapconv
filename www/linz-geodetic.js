
// Module for handling NZ coordinate systems
//
// First cut.  Oddly this ignores (and loses) height information in coordinates...
//
// Namespaces for this module - over complex?

///////////////////////////////////////////////////////////////////////////////////////
// Complex numbers used for NZMG
var Complex = (function () {
    function Complex(re, im) {
        this.re = re;
        this.im = im;
    }
    Complex.prototype.add = function (other) {
        return new Complex(this.re + other.re, this.im + other.im);
    };
    Complex.prototype.subtract = function (other) {
        return new Complex(this.re - other.re, this.im - other.im);
    };
    Complex.prototype.multiply = function (other) {
        return new Complex(this.re * other.re - this.im * other.im, this.re * other.im + this.im * other.re);
    };
    Complex.prototype.divide = function (other) {
        var size = other.re * other.re + other.im * other.im;
        return new Complex((this.re * other.re + this.im * other.im) / size, (this.im * other.re - this.re * other.im) / size);
    };
    Complex.prototype.scale = function (factor) {
        return new Complex(this.re * factor, this.im * factor);
    };
    return Complex;
}())

//////////////////////////////////////////////////////////////////////////
//  NZTM projections
var NZMG = (function () {
    function NZMG() {
        this.a = 6378388.0;
        this.n0 = 6023150.0;
        this.e0 = 2510000.0;
        this.lt0 = -41.0;
        this.ln0 = 173.0;
        this.cfi = [
            0.6399175073,
            -0.1358797613,
            0.063294409,
            -0.02526853,
            0.0117879,
            -0.0055161,
            0.0026906,
            -0.001333,
            0.00067,
            -0.00034
        ];
        this.cfl = [
            1.5627014243,
            0.5185406398,
            -0.03333098,
            -0.1052906,
            -0.0368594,
            0.007317,
            0.01220,
            0.00394,
            -0.0013
        ];
        this.cfb1 = [
            new Complex(0.7557853228, 0.0),
            new Complex(0.249204646, 0.003371507),
            new Complex(-0.001541739, 0.041058560),
            new Complex(-0.10162907, 0.01727609),
            new Complex(-0.26623489, -0.36249218),
            new Complex(-0.6870983, -1.1651967)
        ];
        this.cfb2 = [
            new Complex(1.3231270439, 0.0),
            new Complex(-0.577245789, -0.007809598),
            new Complex(0.508307513, -0.112208952),
            new Complex(-0.15094762, 0.18200602),
            new Complex(1.01418179, 1.64497696),
            new Complex(1.9660549, 2.5127645)
        ];
    }
    NZMG.prototype.toGeodetic = function (coord) {
        var z0 = new Complex((coord[1] - this.n0) / this.a, (coord[0] - this.e0) / this.a);
        var cfb1 = this.cfb1;
        var cfb2 = this.cfb2;
        var z1 = cfb2[5];
        for (var i = 5; i--;)
            z1 = cfb2[i].add(z1.multiply(z0));
        z1 = z1.multiply(z0);
        for (var it = 2; it--;) {
            var zn = cfb1[5].scale(5.0);
            var zd = cfb1[5].scale(6.0);
            for (var i = 4; i; i--) {
                zn = cfb1[i].scale(i).add(zn.multiply(z1));
                zd = cfb1[i].scale(i + 1).add(zd.multiply(z1));
            }
            zn = z0.add(zn.multiply(z1).multiply(z1));
            zd = cfb1[0].add(zd.multiply(z1));
            z1 = zn.divide(zd);
        }
        var ln = this.ln0 + z1.im * 180.0 / Math.PI;
        var cfl = this.cfl;
        var tmp = z1.re;
        var sum = cfl[8];
        for (var i = 8; i--;)
            sum = sum * tmp + cfl[i];
        sum *= tmp / 3600.0e-5;
        var lt = (this.lt0 + sum);
        return [ln, lt];
    };
    NZMG.prototype.toProjection = function (coord) {
        var lt = (coord[1] - this.lt0) * 3600.0e-5;
        var cfi = this.cfi;
        var sum = cfi[9];
        for (var i = 9; i--;)
            sum = sum * lt + cfi[i];
        sum *= lt;
        var cfb1 = this.cfb1;
        var z1 = new Complex(sum, (coord[0] - this.ln0) * Math.PI / 180.0);
        var z0 = cfb1[5];
        for (var i = 5; i--;)
            z0 = cfb1[i].add(z0.multiply(z1));
        z0 = z0.multiply(z1);
        return [this.e0 + z0.im * this.a, this.n0 + z0.re * this.a];
    };
    return NZMG;
}());


//////////////////////////////////////////////////////////////////////////
//  Transverse Mercator projections
var TM =  /** @class */ (function () {
    function TM(a, rf, cm, sf, lto, fe, fn, utom) {
        this.a = a;
        this.rf = rf;
        this.meridian = cm;
        this.scalef = sf;
        this.orglat = lto;
        this.falsee = fe;
        this.falsen = fn;
        this.utom = utom || 1.0;

        var f = rf != 0.0 ? 1.0 / rf : 0.0;
        this.f = f;
        this.e2 = 2.0 * f - f * f;
        this.ep2 = this.e2 / (1.0 - this.e2);
        this.om = this.meridianArc(this.orglat);
    }
    TM.prototype.meridianArc = function (lt) {
        var e2 = this.e2;
        var a = this.a;
        var e4 = e2 * e2;
        var e6 = e4 * e2;

        var A0 = 1 - (e2 / 4.0) - (3.0 * e4 / 64.0) - (5.0 * e6 / 256.0);
        var A2 = (3.0 / 8.0) * (e2 + e4 / 4.0 + 15.0 * e6 / 128.0);
        var A4 = (15.0 / 256.0) * (e4 + 3.0 * e6 / 4.0);
        var A6 = 35.0 * e6 / 3072.0;

        lt = lt * Math.PI / 180.0;
        return a * (A0 * lt - A2 * Math.sin(2 * lt) + A4 * Math.sin(4 * lt) - A6 * Math.sin(6 * lt));
    }
    TM.prototype.footPointLat = function (m) {
        var f = this.f;
        var a = this.a;

        var n = f / (2.0 - f);
        var n2 = n * n;
        var n3 = n2 * n;
        var n4 = n2 * n2;

        var g = a * (1.0 - n) * (1.0 - n2) * (1 + 9.0 * n2 / 4.0 + 225.0 * n4 / 64.0);
        var sig = m / g;

        var phio = sig + (3.0 * n / 2.0 - 27.0 * n3 / 32.0) * Math.sin(2.0 * sig)
            + (21.0 * n2 / 16.0 - 55.0 * n4 / 32.0) * Math.sin(4.0 * sig)
            + (151.0 * n3 / 96.0) * Math.sin(6.0 * sig)
            + (1097.0 * n4 / 512.0) * Math.sin(8.0 * sig);

        return phio;
    }
    TM.prototype.toGeodetic = function (coord) {
        var fn = this.falsen;
        var fe = this.falsee;
        var sf = this.scalef;
        var e2 = this.e2;
        var a = this.a;
        var cm = this.meridian;
        var om = this.om;
        var utom = this.utom;

        var cn1 = (coord[1] - fn) * utom / sf + om;
        var fphi = this.footPointLat(cn1);
        var slt = Math.sin(fphi);
        var clt = Math.cos(fphi);

        var eslt = (1.0 - e2 * slt * slt);
        var eta = a / Math.sqrt(eslt);
        var rho = eta * (1.0 - e2) / eslt;
        var psi = eta / rho;

        var E = (coord[0] - fe) * utom;
        var x = E / (eta * sf);
        var x2 = x * x;

        var t = slt / clt;
        var t2 = t * t;
        var t4 = t2 * t2;

        var trm1 = 1.0 / 2.0;

        var trm2 = ((-4.0 * psi
            + 9.0 * (1 - t2)) * psi
            + 12.0 * t2) / 24.0;

        var trm3 = ((((8.0 * (11.0 - 24.0 * t2) * psi
            - 12.0 * (21.0 - 71.0 * t2)) * psi
            + 15.0 * ((15.0 * t2 - 98.0) * t2 + 15)) * psi
            + 180.0 * ((-3.0 * t2 + 5.0) * t2)) * psi + 360.0 * t4) / 720.0;

        var trm4 = (((1575.0 * t2 + 4095.0) * t2 + 3633.0) * t2 + 1385.0) / 40320.0;

        var lt = fphi + (t * x * E / (sf * rho)) * (((trm4 * x2 - trm3) * x2 + trm2) * x2 - trm1);
        lt *= 180 / Math.PI;

        var trm1 = 1.0;

        var trm2 = (psi + 2.0 * t2) / 6.0;

        var trm3 = (((-4.0 * (1.0 - 6.0 * t2) * psi
            + (9.0 - 68.0 * t2)) * psi
            + 72.0 * t2) * psi
            + 24.0 * t4) / 120.0;

        var trm4 = (((720.0 * t2 + 1320.0) * t2 + 662.0) * t2 + 61.0) / 5040.0;

        var ln = (x / clt) * (((trm4 * x2 - trm3) * x2 + trm2) * x2 - trm1);
        ln = cm - ln * 180.0 / Math.PI;

        return [ln, lt];
    }

    TM.prototype.toProjection = function (coord) {
        var fn = this.falsen;
        var fe = this.falsee;
        var sf = this.scalef;
        var e2 = this.e2;
        var a = this.a;
        var cm = this.meridian;
        var om = this.om;
        var utom = this.utom;

        var ln = coord[0];
        var lt = coord[1];

        var dlon = ln - cm;
        while (dlon > 180.0) dlon -= 360.0;
        while (dlon < -180.0) dlon += 360.0;
        dlon = dlon * Math.PI / 180.0;


        var m = this.meridianArc(lt);
        var slt = Math.sin(lt * Math.PI / 180.0);
        var clt = Math.cos(lt * Math.PI / 180.0);

        var eslt = (1.0 - e2 * slt * slt);
        var eta = a / Math.sqrt(eslt);
        var rho = eta * (1.0 - e2) / eslt;
        var psi = eta / rho;

        var w = dlon;

        var wc = clt * w;
        var wc2 = wc * wc;

        var t = slt / clt;
        var t2 = t * t;
        var t4 = t2 * t2;
        var t6 = t2 * t4;

        var trm1 = (psi - t2) / 6.0;

        var trm2 = (((4.0 * (1.0 - 6.0 * t2) * psi + (1.0 + 8.0 * t2)) * psi - 2.0 * t2) * psi + t4) / 120.0;

        var trm3 = (61 - 479.0 * t2 + 179.0 * t4 - t6) / 5040.0;

        var gce = (sf * eta * dlon * clt) * (((trm3 * wc2 + trm2) * wc2 + trm1) * wc2 + 1.0);

        var ce = gce / utom + fe;

        var trm1 = 1.0 / 2.0;

        var trm2 = ((4.0 * psi + 1) * psi - t2) / 24.0;

        var trm3 = ((((8.0 * (11.0 - 24.0 * t2) * psi
            - 28.0 * (1.0 - 6.0 * t2)) * psi
            + (1.0 - 32.0 * t2)) * psi
            - 2.0 * t2) * psi
            + t4) / 720.0;

        var trm4 = (1385.0 - 3111.0 * t2 + 543.0 * t4 - t6) / 40320.0;

        var gcn = (eta * t) * ((((trm4 * wc2 + trm3) * wc2 + trm2) * wc2 + trm1) * wc2);

        var cn = (gcn + m - om) * sf / utom + fn;

        return [ce, cn];
    }
    return TM;
}());

var NullFunction = (function () {
    function NullFunction() { }

    NullFunction.prototype.toGeodetic = function (x) { return x; },
        NullFunction.prototype.toProjection = function (x) { return x; }
    return NullFunction;
})();



////////////////////////////////////////////////////////////////////////////////
//  Grid based datum conversion using slightly modified NTv2 ascii grid format.
//  Nested grids not supported.
//  Modification are to make the file an exported Javascript string.  Requirement
//  is that each grid data line is the same length.

var GridTransform = (function () {
    function GridTransform(grid_file) {
        this.grid_file = grid_file;
        this.loaded = false;
        this.load_failed = false;

    }

    GridTransform.prototype.load = function () {
        if (this.loaded) return true;
        if (this.load_failed) return false;
        var gridname = this.grid_file;
        gridname = gridname.replace(/^.*[\\\/]/, '');
        gridname = gridname.replace(/\..*/, '');
        this.load_failed = true;
        try {
            griddata = Geodetic.GridData[gridname];
        }
        catch (err) {
        }
        if (!griddata) return false;
        this.lon0 = griddata.lon0
        this.lat0 = griddata.lat0;
        this.lon1 = griddata.lon0;
        this.lat1 = griddata.lat0;
        this.loninc = griddata.loninc;
        this.latinc = griddata.latinc;
        this.nlat = griddata.nlat;
        this.nlon = griddata.nlon;
        this.datascale = griddata.datascale;
        nlon = this.nlon;
        nlat = this.nlat;
        this.getGridData = function (ilon, ilat) {
            if (ilon < 0 || ilon >= this.nlon) return;
            if (ilat < 0 || ilat >= this.nlat) return;
            return griddata.data[ilon][ilat];
        }
        this.loaded = true;
        this.load_failed = false;
        return this.loaded;

    }

    GridTransform.prototype.calcOffset = function (coord) {
        var glon = (coord[0] - this.lon0) / this.loninc;
        var glat = (coord[1] - this.lat0) / this.latinc;
        var ilon = Math.floor(glon);
        var ilat = Math.floor(glat);
        if (ilon < 0 || ilon > this.nlon - 2) return;
        if (ilat < 0 || ilat > this.nlat - 2) return;
        var d00 = this.getGridData(ilon, ilat)
        var d10 = this.getGridData(ilon + 1, ilat)
        var d01 = this.getGridData(ilon, ilat + 1)
        var d11 = this.getGridData(ilon + 1, ilat + 1)
        glon -= ilon;
        glat -= ilat;
        offset = [
            ((d00[0] * (1 - glon) + d10[0] * glon) * (1 - glat) + (d01[0] * (1 - glon) + d11[0] * glon) * glat) * this.datascale,
            ((d00[1] * (1 - glon) + d10[1] * glon) * (1 - glat) + (d01[1] * (1 - glon) + d11[1] * glon) * glat) * this.datascale
        ];
        return offset;
    }

    GridTransform.prototype.toGeodetic = function (coord) {
        if (!this.load()) return;
        var offset = this.calcOffset(coord);
        if (!offset) return;
        return [coord[0] + offset[0], coord[1] + offset[1]];
    }

    GridTransform.prototype.toProjection = function (coord) {
        if (!this.load()) return;
        var offset = this.calcOffset(coord);
        if (!offset) return;
        var coord2 = [coord[0] - offset[0], coord[1] - offset[1]];
        offset = this.calcOffset(coord2);
        return [coord[0] - offset[0], coord[1] - offset[1]];
    }
    return GridTransform;
})();

//////////////////////////////////////////////////////////////////////////////////////////////
//  Coordinate systems - hard-coded for the moment!
var CoordSys = /** @class */ (function () {
    function CoordSys(code, name, type, baseCode, transform) {
        this.code = code;
        this.name = name;
        this.type = type;
        this.baseCode = baseCode;
        this.transform = transform;
    }
    CoordSys.getCoordSys = function (code) {
        var GEOCENTRIC = 0; // XYZ coords
        var GEODETIC = 1; // Lat/lon/ellipsoidal height
        var PROJECTION = 2; // Easting/northing
        var CoordSysList = {
            'NZGD2000': new CoordSys('NZGD2000', 'New Zealand Geodetic Datum 2000', GEODETIC, null, null),
            'NZTM': new CoordSys('NZTM', 'New Zealand Transverse Mercator', PROJECTION, 'NZGD2000', new TM(6378137.0, 298.257222, 173.0, 0.9996, 0.0, 1600000.0, 10000000.0, 1.0)),
            'NZGD1949': new CoordSys('NZGD1949', 'New Zealand Geodetic Datum 1949', GEODETIC, 'NZGD2000', new GridTransform('./nzgd2kgrid9911.asc.json')),
            'NZMG': new CoordSys('NZMG', 'New Zealand Map Grid', PROJECTION, 'NZGD1949', new NZMG())
        };
        return CoordSysList[code];
    };
    CoordSys.MaxDepth = 10;
    return CoordSys;
}());


//////////////////////////////////////////////////////////////////////////
// Location: Represents a physical location, which may be expressed in a 
// number of coordinate systems.  It is defined with a specified coordinate
// system.  When another is requested it works out a transformation path to 
// the new system.  Each intermediate coordinate, as well as the final 
// coordinate, is saved with the location
var Location = (function () {

    function Location(coordSysCode, coordinates) {
        this.coordinates = {};
        this.coordinates[coordSysCode] = coordinates;
    }

    // Express the location in a selected coordinate system.
    Location.prototype.as = function (coordSysCode) {
        // If already available, just return it.

        if (coordSysCode in this.coordinates) return this.coordinates[coordSysCode];
        var getcrdsys = CoordSys.getCoordSys;
        var coordsys = getcrdsys(coordSysCode);
        if (!coordsys) return;

        // Find options for transforming to required system
        var crdSysStages = [coordsys.code];
        var transformStages = [[coordsys.transform, 0]];
        var gotTransform = false;
        var depth = CoordSys.MaxDepth;
        while (coordsys && depth-- > 0) {
            if (coordsys.baseCode in this.coordinates) {
                crdSysStages.unshift(coordsys.baseCode);
                transformStages.unshift(null)
                gotTransform = true;
                break;
            }
            var coordsys = getcrdsys(coordsys.baseCode);
            if (!coordsys) break;
            crdSysStages.unshift(coordsys.code);
            transformStages.unshift([coordsys.transform, 0]);
        }

        if (!gotTransform) {
            var targetCrdSys = crdSysStages;
            var targetTransforms = transformStages;

            // Find candidate to transform from
            for (cscode in this.coordinates) {
                fcs = getcrdsys(cscode);
                // Only interested in coordinate systems with base codes...
                if (!fcs || !fcs.baseCode || fcs.baseCode in this.coordinates) continue;
                crdSysStages = [cscode, fcs.baseCode];
                transformStages = [null, [fcs.transform, 1]];
                depth = CoordSys.MaxDepth;

                while (depth-- > 0) {
                    targetIndex = targetCrdSys.indexOf(fcs.baseCode);
                    if (targetIndex >= 0) {
                        crdSysStages = crdSysStages.concat(targetCrdSys.slice(targetIndex + 1));
                        transformStages = transformStages.concat(targetTransforms.slice(targetIndex + 1));
                        gotTransform = true;
                        break;
                    }
                    fcs = getcrdsys(fcs.baseCode);
                    if (!fcs || !fcs.baseCode || fcs.baseCode in this.coordinates) break;
                    crdSysStages.push(fcs.baseCode);
                    transformStages.push([fcs.transform, 1]);
                }
                if (gotTransform) break;
            }
        }
        if (!gotTransform) return;

        var coord = this.coordinates[crdSysStages[0]];
        for (i = 1; i < crdSysStages.length; i++) {
            tfm = transformStages[i];
            coord = tfm[1] ? tfm[0].toGeodetic(coord) : tfm[0].toProjection(coord);
            this.coordinates[crdSysStages[i]] = coord;
        }
        return this.coordinates[coordSysCode];
    }
    return Location;
})();