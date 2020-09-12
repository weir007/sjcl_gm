/**
 * base class for all ecc operations.
 * @namespace
 */
sjcl.ecc = {};

/**
 * Represents a point on a curve in affine coordinates.
 * @constructor
 * @param {sjcl.ecc.curve} curve The curve that this point lies on.
 * @param {bigInt} x The x coordinate.
 * @param {bigInt} y The y coordinate.
 */
sjcl.ecc.point = function(curve,x,y) {
  if (x === undefined) {
    this.isIdentity = true;
  } else {
    if (x instanceof sjcl.bn) {
      x = new curve.field(x);
    }
    if (y instanceof sjcl.bn) {
      y = new curve.field(y);
    }

    this.x = x;
    this.y = y;

    this.isIdentity = false;
  }
  this.curve = curve;
};



sjcl.ecc.point.prototype = {
  toJac: function() {
    return new sjcl.ecc.pointJac(this.curve, this.x, this.y, new this.curve.field(1));
  },

  multG: function(k) {
	return this.toJac().multG(k).toAffine();
	//return this.toJac().mult(k, this).toAffine();//close optimization
  },
  
  mult: function(k) {
    return this.toJac().mult(k, this).toAffine();
  },

  /**
   * Multiply this point by k, added to affine2*k2, and return the answer in Jacobian coordinates.
   * @param {bigInt} k The coefficient to multiply this by.
   * @param {bigInt} k2 The coefficient to multiply affine2 this by.
   * @param {sjcl.ecc.point} affine The other point in affine coordinates.
   * @return {sjcl.ecc.pointJac} The result of the multiplication and addition, in Jacobian coordinates.
   */
  mult2: function(k, k2, affine2) {
    return this.toJac().mult2(k, this, k2, affine2).toAffine();
  },
  /* precompute window=4 */
  multiples: function() {
    var m, i, j;
    if (this._multiples === undefined) {
      j = this.toJac().doubl();
      m = this._multiples = [new sjcl.ecc.point(this.curve), this, j.toAffine()];
      for (i=3; i<16; i++) {
        j = j.add(this);
        m.push(j.toAffine());
      }
    }
    return this._multiples;
  },

  negate: function() {
    var newY = new this.curve.field(0).sub(this.y).normalize().reduce();
    return new sjcl.ecc.point(this.curve, this.x, newY);
  },

  isValid: function() {
    return this.y.square().equals(this.curve.b.add(this.x.mul(this.curve.a.add(this.x.square()))));
  },

  toBits: function() {
    return sjcl.bitArray.concat(this.x.toBits(), this.y.toBits());
  }
};

/**
 * Represents a point on a curve in Jacobian coordinates. Coordinates can be specified as bigInts or strings (which
 * will be converted to bigInts).
 *
 * @constructor
 * @param {bigInt/string} x The x coordinate.
 * @param {bigInt/string} y The y coordinate.
 * @param {bigInt/string} z The z coordinate.
 * @param {sjcl.ecc.curve} curve The curve that this point lies on.
 */
sjcl.ecc.pointJac = function(curve, x, y, z) {
  if (x === undefined) {
    this.isIdentity = true;
  } else {
    this.x = x;
    this.y = y;
    this.z = z;
    this.isIdentity = false;
  }
  this.curve = curve;
};

sjcl.ecc.pointJac.prototype = {
  /**
   * Adds S and T and returns the result in Jacobian coordinates. Note that S must be in Jacobian coordinates and T must be in affine coordinates.
   * @param {sjcl.ecc.pointJac} S One of the points to add, in Jacobian coordinates.
   * @param {sjcl.ecc.point} T The other point to add, in affine coordinates.
   * @return {sjcl.ecc.pointJac} The sum of the two points, in Jacobian coordinates.
   */
  add: function(T) {
    var S = this, sz2, c, d, c2, x1, x2, x, y1, y2, y, z;
    if (S.curve !== T.curve) {
      throw new sjcl.exception.invalid("sjcl.ecc.add(): Points must be on the same curve to add them!");
    }

    if (S.isIdentity) {
      return T.toJac();
    } else if (T.isIdentity) {
      return S;
    }

    sz2 = S.z.square();
    c = T.x.mul(sz2).subM(S.x);

    if (c.equals(0)) {
      if (S.y.equals(T.y.mul(sz2.mul(S.z)))) {
        // same point
        return S.doubl();
      } else {
        // inverses
        return new sjcl.ecc.pointJac(S.curve);
      }
    }

    d = T.y.mul(sz2.mul(S.z)).subM(S.y);
    c2 = c.square();

    x1 = d.square();
    x2 = c.square().mul(c).addM( S.x.add(S.x).mul(c2) );
    x  = x1.subM(x2);

    y1 = S.x.mul(c2).subM(x).mul(d);
    y2 = S.y.mul(c.square().mul(c));
    y  = y1.subM(y2);

    z  = S.z.mul(c);

    return new sjcl.ecc.pointJac(this.curve,x,y,z);
  },

  /**
   * doubles this point.
   * @return {sjcl.ecc.pointJac} The doubled point.
   */
  doubl: function() {
    if (this.isIdentity) { return this; }

    var
      y2 = this.y.square(),
      a  = y2.mul(this.x.mul(4)),
      b  = y2.square().mul(8),
      z2 = this.z.square(),
      c  = this.curve.a.toString() == (new sjcl.bn(-3)).toString() ?
                this.x.sub(z2).mul(3).mul(this.x.add(z2)) :
                this.x.square().mul(3).add(z2.square().mul(this.curve.a)),
      x  = c.square().subM(a).subM(a),
      y  = a.sub(x).mul(c).subM(b),
      z  = this.y.add(this.y).mul(this.z);
    return new sjcl.ecc.pointJac(this.curve, x, y, z);
  },

  /**
   * Returns a copy of this point converted to affine coordinates.
   * @return {sjcl.ecc.point} The converted point.
   */
  toAffine: function() {
    if (this.isIdentity || this.z.equals(0)) {
      return new sjcl.ecc.point(this.curve);
    }
    var zi = this.z.inverse(), zi2 = zi.square();
    return new sjcl.ecc.point(this.curve, this.x.mul(zi2).fullReduce(), this.y.mul(zi2.mul(zi)).fullReduce());
  },

  /* window = 4 */
  multG: function(k) {
	if (typeof(k) === "number") {
      k = [k];
    } else if (k.limbs !== undefined) {
      k = k.normalize().limbs;
    }
	while (k.length < this.curve.r.limbs.length) k.push(0);
	var i, j, t, r = sjcl.bn.prototype.radix, out = new sjcl.ecc.point(this.curve).toJac();
	for (i=((this.curve.r.bitLength()+3)>>2)-1; i>=0; i--) {
	  t  = (k[Math.floor(i/r)] >> (i%r)) & 1;
	  t += ((k[Math.floor((i+64)/r)]  >> ((i+64)%r))  & 1) << 1;
	  t += ((k[Math.floor((i+128)/r)] >> ((i+128)%r)) & 1) << 2;
	  t += ((k[Math.floor((i+192)/r)] >> ((i+192)%r)) & 1) << 3;
	  out = out.doubl().add(this.curve.mG[t]);
    }
	return out;
  },
  
  /**
   * Multiply this point by k and return the answer in Jacobian coordinates.
   * @param {bigInt} k The coefficient to multiply by.
   * @param {sjcl.ecc.point} affine This point in affine coordinates.
   * @return {sjcl.ecc.pointJac} The result of the multiplication, in Jacobian coordinates.
   */
  mult: function(k, affine) {
    if (typeof(k) === "number") {
      k = [k];
    } else if (k.limbs !== undefined) {
      k = k.normalize().limbs;
    }

    var i, j, t, out = new sjcl.ecc.point(this.curve).toJac(), multiples = affine.multiples();

	/*
	//it's not constant-time for k not filled to fix length (this.curve.r.bitLength() -> this.curve.r.limbs.length)
    for (i=k.length-1; i>=0; i--) {
      for (j=sjcl.bn.prototype.radix-4; j>=0; j-=4) {
        out = out.doubl().doubl().doubl().doubl().add(multiples[k[i]>>j & 0xF]);
      }
    } */
	
	/* Const-time scalar-mult */
	for (i=this.curve.r.limbs.length-1; i>=0; i--) {
	  t = i < k.length ? k[i] : 0;
      for (j=sjcl.bn.prototype.radix-4; j>=0; j-=4) {
        out = out.doubl().doubl().doubl().doubl().add(multiples[t>>j & 0xF]);
      }
    }
	return out;
  },  
  /*
  mult: function(k, affine, st=true) {
    if (typeof(k) === "number") {
      k = [k];
    } else if (k.limbs !== undefined) {
      k = k.normalize().limbs;
    }

    var i, j, t, out = new sjcl.ecc.point(this.curve).toJac(), multiples = affine.multiples();

	/* 
	//it's not constant-time for k not filled to fix length (this.curve.r.bitLength() -> this.curve.r.limbs.length)
    for (i=k.length-1; i>=0; i--) {
      for (j=sjcl.bn.prototype.radix-4; j>=0; j-=4) {
        out = out.doubl().doubl().doubl().doubl().add(multiples[k[i]>>j & 0xF]);
      }
    } * /
	
	/* Const-time scalar-mult * /
	if(typeof Worker === "undefined" || typeof window === "undefined" || st) {
      for (i=this.curve.r.limbs.length-1; i>=0; i--) {
	    t = i < k.length ? k[i] : 0;
        for (j=sjcl.bn.prototype.radix-4; j>=0; j-=4) {
          out = out.doubl().doubl().doubl().doubl().add(multiples[t>>j & 0xF]);
        }
      }
	}
	else{
	  //var mul_s = function (res, muls = multiples, w = Math.ceil(this.curve.r.length/x), k, s, e = s+w, r = sjcl.bn.prototype.radix) {
	  var mul_s = function mul_s () {
	    self.addEventListener('message', function(evt) {
	  	/* receive params * /
	  	var k = evt.data['k'],
      		muls = evt.data['muls'],
	  		res = evt.data['res'],
	  		s = evt.data['s'],
	  		e = evt.data['e'],
			r = evt.data['r'],
	  	    i, j, t;
	       for (i=s; i<e; i++) {
	         t = typeof k[i] === 'undefined' ? 0 : k[i];
	         for (j=r; j>=0; j-=4) {
	       	  res = res.doubl().doubl().doubl().doubl().add(muls[t>>j & 0xF]);
	         }
	       }
           postMessage(res);
	    })
      };
	  
	  var x = 2, scalar = new sjcl.bn(2), res = [out], blob, blobURL, worker = [], w = Math.ceil(this.curve.r.limbs.length / x);
	  
	  for (i=0; i<x; i++) {
		if (i>0) {
		  scalar = scalar.power(i*w*sjcl.bn.prototype.radix);
		  res[i] = this.mult(scalar, affine, true);
	  	}
		blob = new Blob([mul_s+"\nmul_s()"]);
        blobURL = window.URL.createObjectURL(blob);
        worker[i] = new Worker(blobURL);
	    worker[i].onmessage = function(e) {
		  //console.log('th_'+i+' : '+e.data);
		  out = out.add(e.data); 
		};
		//var postMessageTemp = window.postMessage;
        //window.postMessage = function(message, targetOrigin, transfer)
        //{
          //postMessageTemp(JSON.parse(JSON.stringify(message)), targetOrigin, transfer)
        //};
	    //worker[i].postMessage({'res':res[i], 'k':k, 'muls':multiples, 's':i*w, 'e':i*w+w});
	    
		worker[i].postMessage({'res':res[i], 'k':k, 's':i*w, 'e':i*w+w, 'r':sjcl.bn.prototype.radix},[res[i]]);
	  }
	}
    return out;
  },
*/
  /**
   * Multiply this point by k, added to affine2*k2, and return the answer in Jacobian coordinates.
   * @param {bigInt} k The coefficient to multiply this by.
   * @param {sjcl.ecc.point} affine This point in affine coordinates.
   * @param {bigInt} k2 The coefficient to multiply affine2 this by.
   * @param {sjcl.ecc.point} affine The other point in affine coordinates.
   * @return {sjcl.ecc.pointJac} The result of the multiplication and addition, in Jacobian coordinates.
   */
  mult2: function(k1, affine, k2, affine2) {
    if (typeof(k1) === "number") {
      k1 = [k1];
    } else if (k1.limbs !== undefined) {
      k1 = k1.normalize().limbs;
    }

    if (typeof(k2) === "number") {
      k2 = [k2];
    } else if (k2.limbs !== undefined) {
      k2 = k2.normalize().limbs;
    }

    var i, j, out = new sjcl.ecc.point(this.curve).toJac(), m1 = affine.multiples(),
        m2 = affine2.multiples(), l1, l2;

    for (i=Math.max(k1.length,k2.length)-1; i>=0; i--) {
      l1 = k1[i] | 0;
      l2 = k2[i] | 0;
      for (j=sjcl.bn.prototype.radix-4; j>=0; j-=4) {
        out = out.doubl().doubl().doubl().doubl().add(m1[l1>>j & 0xF]).add(m2[l2>>j & 0xF]);
      }
    }

    return out;
  },

  negate: function() {
    return this.toAffine().negate().toJac();
  },

  isValid: function() {
    var z2 = this.z.square(), z4 = z2.square(), z6 = z4.mul(z2);
    return this.y.square().equals(
             this.curve.b.mul(z6).add(this.x.mul(
               this.curve.a.mul(z4).add(this.x.square()))));
  }
};

/**
 * Construct an elliptic curve. Most users will not use this and instead start with one of the NIST curves defined below.
 *
 * @constructor
 * @param {bigInt} p The prime modulus.
 * @param {bigInt} r The prime order of the curve.
 * @param {bigInt} a The constant a in the equation of the curve y^2 = x^3 + ax + b (for NIST curves, a is always -3).
 * @param {bigInt} x The x coordinate of a base point of the curve.
 * @param {bigInt} y The y coordinate of a base point of the curve.
 */
sjcl.ecc.curve = function(Field, r, a, b, x, y) {
  this.field = Field;
  this.r = new sjcl.bn(r);
  this.a = new Field(a);
  this.b = new Field(b);
  this.G = new sjcl.ecc.point(this, new Field(x), new Field(y));
  /* maybe I can precompute something here */
  var m = this.mG = [new sjcl.ecc.point(this), this.G];
  var s = new sjcl.bn(2).power((this.r.bitLength()+3)>>2);
  m[2] = m[1].mult(s);
  m[4] = m[2].mult(s);
  m[8] = m[4].mult(s);
  
  for (var i=3; i<16; i++)
    m[i] = m[i&1].toJac().add(m[i&2]).add(m[i&4]).add(m[i&8]).toAffine();
};

sjcl.ecc.curve.prototype.fromBits = function (bits) {
  var w = sjcl.bitArray, l = this.field.prototype.exponent + 7 & -8,
      p = new sjcl.ecc.point(this, this.field.fromBits(w.bitSlice(bits, 0, l)),
                             this.field.fromBits(w.bitSlice(bits, l, 2*l)));
  if (!p.isValid()) {
    throw new sjcl.exception.corrupt("not on the curve!");
  }
  return p;
};

sjcl.ecc.curves = {
  c192: new sjcl.ecc.curve(
    sjcl.bn.prime.p192,
    "0xffffffffffffffffffffffff99def836146bc9b1b4d22831",
    -3,
    "0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1",
    "0x188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012",
    "0x07192b95ffc8da78631011ed6b24cdd573f977a11e794811"),

  c224: new sjcl.ecc.curve(
    sjcl.bn.prime.p224,
    "0xffffffffffffffffffffffffffff16a2e0b8f03e13dd29455c5c2a3d",
    -3,
    "0xb4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4",
    "0xb70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21",
    "0xbd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34"),

  c256: new sjcl.ecc.curve(
    sjcl.bn.prime.p256,
    "0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551",
    -3,
    "0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
    "0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296",
    "0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5"),

  c384: new sjcl.ecc.curve(
    sjcl.bn.prime.p384,
    "0xffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973",
    -3,
    "0xb3312fa7e23ee7e4988e056be3f82d19181d9c6efe8141120314088f5013875ac656398d8a2ed19d2a85c8edd3ec2aef",
    "0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7",
    "0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f"),
    
  c521: new sjcl.ecc.curve(
    sjcl.bn.prime.p521,
    "0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409",
    -3,
    "0x051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00",
    "0xC6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66",
    "0x11839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650"),

  k192: new sjcl.ecc.curve(
    sjcl.bn.prime.p192k,
    "0xfffffffffffffffffffffffe26f2fc170f69466a74defd8d",
    0,
    3,
    "0xdb4ff10ec057e9ae26b07d0280b7f4341da5d1b1eae06c7d",
    "0x9b2f2f6d9c5628a7844163d015be86344082aa88d95e2f9d"),

  k224: new sjcl.ecc.curve(
    sjcl.bn.prime.p224k,
    "0x010000000000000000000000000001dce8d2ec6184caf0a971769fb1f7",
    0,
    5,
    "0xa1455b334df099df30fc28a169a467e9e47075a90f7e650eb6b7a45c",
    "0x7e089fed7fba344282cafbd6f7e319f7c0b0bd59e2ca4bdb556d61a5"),

  k256: new sjcl.ecc.curve(
    sjcl.bn.prime.p256k,
    "0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",
    0,
    7,
    "0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
    "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8"),

  sm2: new sjcl.ecc.curve(/* sm2 */
    sjcl.bn.prime.p256s,
    "0xfffffffeffffffffffffffffffffffff7203df6b21c6052b53bbf40939d54123",
	"0xfffffffeffffffffffffffffffffffffffffffff00000000fffffffffffffffc",
	"0x28e9fa9e9d9f5e344d5a9e4bcf6509a7f39789f515ab8f92ddbcbd414d940e93",
	"0x32c4ae2c1f1981195f9904466a39c9948fe30bbff2660be1715a4589334c74c7",
	"0xbc3736a2f4f6779c59bdcee36b692153d0a9877cc62a474002df32e52139f0a0")
};

sjcl.ecc.curveName = function (curve) {
  var curcurve;
  for (curcurve in sjcl.ecc.curves) {
    if (sjcl.ecc.curves.hasOwnProperty(curcurve)) {
      if (sjcl.ecc.curves[curcurve] === curve) {
        return curcurve;
      }
    }
  }

  throw new sjcl.exception.invalid("no such curve");
};

sjcl.ecc.deserialize = function (key) {
  var types = ["elGamal", "ecdsa", "sm2"];//weir007

  if (!key || !key.curve || !sjcl.ecc.curves[key.curve]) { throw new sjcl.exception.invalid("invalid serialization"); }
  if (types.indexOf(key.type) === -1) { throw new sjcl.exception.invalid("invalid type"); }

  var curve = sjcl.ecc.curves[key.curve];

  if (key.secretKey) {
    if (!key.exponent) { throw new sjcl.exception.invalid("invalid exponent"); }
    var exponent = new sjcl.bn(key.exponent);
    return new sjcl.ecc[key.type].secretKey(curve, exponent);
  } else {
    if (!key.point) { throw new sjcl.exception.invalid("invalid point"); }
    
    var point = curve.fromBits(sjcl.codec.hex.toBits(key.point));
    return new sjcl.ecc[key.type].publicKey(curve, point);
  }
};

/** our basicKey classes
*/
sjcl.ecc.basicKey = {
  /** ecc publicKey.
   * @constructor
   * @param {curve} curve the elliptic curve
   * @param {point} point the point on the curve
   */
  publicKey: function(curve, point) {
    this._curve = curve;
    this._curveBitLength = curve.r.bitLength();
    if (point instanceof Array) {
      this._point = curve.fromBits(point);
    } else {
      this._point = point;
    }

    this.serialize = function () {
      var curveName = sjcl.ecc.curveName(curve);
      return {
        type: this.getType(),
        secretKey: false,
        point: sjcl.codec.hex.fromBits(this._point.toBits()),
        curve: curveName
      };
    };

    /** get this keys point data
    * @return x and y as bitArrays
    */
    this.get = function() {
      var pointbits = this._point.toBits();
      var len = sjcl.bitArray.bitLength(pointbits);
      var x = sjcl.bitArray.bitSlice(pointbits, 0, len/2);
      var y = sjcl.bitArray.bitSlice(pointbits, len/2);
      return { x: x, y: y };
    };
  },

  /** ecc secretKey
  * @constructor
  * @param {curve} curve the elliptic curve
  * @param exponent
  */
  secretKey: function(curve, exponent) {
    this._curve = curve;
    this._curveBitLength = curve.r.bitLength();
    this._exponent = exponent;

    this.serialize = function () {
      var exponent = this.get();
      var curveName = sjcl.ecc.curveName(curve);
      return {
        type: this.getType(),
        secretKey: true,
        exponent: sjcl.codec.hex.fromBits(exponent),
        curve: curveName
      };
    };

    /** get this keys exponent data
    * @return {bitArray} exponent
    */
    this.get = function () {
      return this._exponent.toBits();
    };
  }
};

/** @private * /
sjcl.ecc.basicKey.generateKeys = function(cn) {
  return function generateKeys(curve, paranoia, sec) {
    curve = curve || 256;

    if (typeof curve === "number") {
      curve = sjcl.ecc.curves['c'+curve];
      if (curve === undefined) {
        throw new sjcl.exception.invalid("no such curve");
      }
    }
    sec = sec || sjcl.bn.random(curve.r, paranoia);

    var pub = curve.G.mult(sec);
    return { pub: new sjcl.ecc[cn].publicKey(curve, pub),
             sec: new sjcl.ecc[cn].secretKey(curve, sec) };
  };
};
*/

sjcl.ecc.basicKey.generateKeys = function(cn) {
  return function generateKeys(curve, paranoia, sec) {
    if (cn === 'sm2')
	{
	  curve = sjcl.ecc.curves[cn];
	}
	else
	{
      curve |= 256;
	  if (typeof curve === "number") {
        curve = sjcl.ecc.curves['c'+curve];
        if (curve === undefined) {
          throw new sjcl.exception.invalid("no such curve");
        }
      }
	}
	sec = sec || sjcl.bn.random(curve.r, paranoia);
    //var pub = curve.G.mult(sec);//multG
    var pub = curve.G.multG(sec);
    return { pub: new sjcl.ecc[cn].publicKey(curve, pub),
             sec: new sjcl.ecc[cn].secretKey(curve, sec) };
  };
};

/** elGamal keys */
sjcl.ecc.elGamal = {
  /** generate keys
  * @function
  * @param curve
  * @param {int} paranoia Paranoia for generation (default 6)
  * @param {secretKey} sec secret Key to use. used to get the publicKey for ones secretKey
  */
  generateKeys: sjcl.ecc.basicKey.generateKeys("elGamal"),
  /** elGamal publicKey.
  * @constructor
  * @augments sjcl.ecc.basicKey.publicKey
  */
  publicKey: function (curve, point) {
    sjcl.ecc.basicKey.publicKey.apply(this, arguments);
  },
  /** elGamal secretKey
  * @constructor
  * @augments sjcl.ecc.basicKey.secretKey
  */
  secretKey: function (curve, exponent) {
    sjcl.ecc.basicKey.secretKey.apply(this, arguments);
  }
};

sjcl.ecc.elGamal.publicKey.prototype = {
  /** Kem function of elGamal Public Key
  * @param paranoia paranoia to use for randomization.
  * @return {object} key and tag. unkem(tag) with the corresponding secret key results in the key returned.
  */
  kem: function(paranoia) {
    var sec = sjcl.bn.random(this._curve.r, paranoia),
        tag = this._curve.G.mult(sec).toBits(),
        key = sjcl.hash.sha256.hash(this._point.mult(sec).toBits());
    return { key: key, tag: tag };
  },
  
  getType: function() {
    return "elGamal";
  }
};

sjcl.ecc.elGamal.secretKey.prototype = {
  /** UnKem function of elGamal Secret Key
  * @param {bitArray} tag The Tag to decrypt.
  * @return {bitArray} decrypted key.
  */
  unkem: function(tag) {
    return sjcl.hash.sha256.hash(this._curve.fromBits(tag).mult(this._exponent).toBits());
  },

  /** Diffie-Hellmann function
  * @param {elGamal.publicKey} pk The Public Key to do Diffie-Hellmann with
  * @return {bitArray} diffie-hellmann result for this key combination.
  */
  dh: function(pk) {
    return sjcl.hash.sha256.hash(pk._point.mult(this._exponent).toBits());
  },

  /** Diffie-Hellmann function, compatible with Java generateSecret
  * @param {elGamal.publicKey} pk The Public Key to do Diffie-Hellmann with
  * @return {bitArray} undigested X value, diffie-hellmann result for this key combination,
  * compatible with Java generateSecret().
  */
  dhJavaEc: function(pk) {
    return pk._point.mult(this._exponent).x.toBits();
  }, 

  getType: function() {
    return "elGamal";
  }
};

/** ecdsa keys */
sjcl.ecc.ecdsa = {
  /** generate keys
  * @function
  * @param curve
  * @param {int} paranoia Paranoia for generation (default 6)
  * @param {secretKey} sec secret Key to use. used to get the publicKey for ones secretKey
  */
  generateKeys: sjcl.ecc.basicKey.generateKeys("ecdsa")
};

/** ecdsa publicKey.
* @constructor
* @augments sjcl.ecc.basicKey.publicKey
*/
sjcl.ecc.ecdsa.publicKey = function (curve, point) {
  sjcl.ecc.basicKey.publicKey.apply(this, arguments);
};

/** specific functions for ecdsa publicKey. */
sjcl.ecc.ecdsa.publicKey.prototype = {
  /** Diffie-Hellmann function
  * @param {bitArray} hash hash to verify.
  * @param {bitArray} rs signature bitArray.
  * @param {boolean}  fakeLegacyVersion use old legacy version
  */
  verify: function(hash, rs, fakeLegacyVersion) {
    if (sjcl.bitArray.bitLength(hash) > this._curveBitLength) {
      hash = sjcl.bitArray.clamp(hash, this._curveBitLength);
    }
    var w = sjcl.bitArray,
        R = this._curve.r,
        l = this._curveBitLength,
        r = sjcl.bn.fromBits(w.bitSlice(rs,0,l)),
        ss = sjcl.bn.fromBits(w.bitSlice(rs,l,2*l)),
        s = fakeLegacyVersion ? ss : ss.inverseMod(R),
        hG = sjcl.bn.fromBits(hash).mul(s).mod(R),
        hA = r.mul(s).mod(R),
        r2 = this._curve.G.mult2(hG, hA, this._point).x;
    if (r.equals(0) || ss.equals(0) || r.greaterEquals(R) || ss.greaterEquals(R) || !r2.equals(r)) {
      if (fakeLegacyVersion === undefined) {
        return this.verify(hash, rs, true);
      } else {
        throw (new sjcl.exception.corrupt("signature didn't check out"));
      }
    }
    return true;
  },

  getType: function() {
    return "ecdsa";
  }
};

/** ecdsa secretKey
* @constructor
* @augments sjcl.ecc.basicKey.publicKey
*/
sjcl.ecc.ecdsa.secretKey = function (curve, exponent) {
  sjcl.ecc.basicKey.secretKey.apply(this, arguments);
};

/** specific functions for ecdsa secretKey. */
sjcl.ecc.ecdsa.secretKey.prototype = {
  /** Diffie-Hellmann function
  * @param {bitArray} hash hash to sign.
  * @param {int} paranoia paranoia for random number generation
  * @param {boolean} fakeLegacyVersion use old legacy version
  */
  sign: function(hash, paranoia, fakeLegacyVersion, fixedKForTesting) {
    if (sjcl.bitArray.bitLength(hash) > this._curveBitLength) {
      hash = sjcl.bitArray.clamp(hash, this._curveBitLength);
    }
    var R  = this._curve.r,
        l  = R.bitLength(),
        k  = fixedKForTesting || sjcl.bn.random(R.sub(1), paranoia).add(1),
        //r  = this._curve.G.mult(k).x.mod(R),
        r  = this._curve.G.multG(k).x.mod(R),
        ss = sjcl.bn.fromBits(hash).add(r.mul(this._exponent)),
        s  = fakeLegacyVersion ? ss.inverseMod(R).mul(k).mod(R)
             : ss.mul(k.inverseMod(R)).mod(R);
    return sjcl.bitArray.concat(r.toBits(l), s.toBits(l));
  },

  getType: function() {
    return "ecdsa";
  }
};


/* test */
function dump(a, tag) {
  var str = tag+':';
  for (var i=0; i<a.length; i++)
    str = str.concat((a[i]<0 ? a[i]+0x100000000 : a[i]).toString(16)+' ');
  console.log(str);
}

function dumpb(a, tag) {
  var str = tag+':';
  for (var i=0; i<a.length; i++)
    str = str.concat((a[i]<0 ? a[i]+0x100 : a[i]).toString(16)+' ');
  console.log(str);
}

/**
 * SM2
 * Including Key Generation, Sign/Verify and Encryption/Decryption,
 * consequently, KDF based on sm3 is needed.
 */
/** sm2 keys */
sjcl.ecc.sm2 = {
  /** generate keys
  * @function
  * @param curve
  * @param {int} paranoia Paranoia for generation (default 6)
  * @param {secretKey} sec secret Key to use. used to get the publicKey for ones secretKey
  */
  //generateKeys: sjcl.ecc.basicKey.generateKeys("sm2"),
  generateKeys: function () {
	var key = sjcl.ecc.basicKey.generateKeys("sm2")();
	/* set pubkey, for signing */
	//Object.assign(key.sec, {'_pk':key.pub._point});
	key.sec._pk = key.pub._point;
	return key;
  },
  /**
   * ZA = H(ENTLA || IDA || a || b || xG || yG || xA || yA)
   * @param {ecc.curve} curve 
   * @param {xcoordinate} px of public key
   * @param {ycoordinate} py of public key
   * @param {string | bitArray} ida UserID
   * @param {sjcl.hash.xxx.hash} hash hash function   
   */
  HM : function (curve, px, py, ida, hash) {
    ida = ida || "1234567812345678";
	hash = hash || sjcl.hash.sm3.hash;
	
	var bytes = sjcl.codec.bytes,
	  string = sjcl.codec.utf8String,
	  hex = sjcl.codec.hex;
	if (typeof ida === 'string')
	  ida = bytes.fromBits(string.toBits(ida));
  
	var x = bytes.fromBits(hex.toBits(px.toString(16))), y = bytes.fromBits(hex.toBits(py.toString(16))),
        la = ida.length << 3, 
	    entla = [(la>>8) & 0xff, la & 0xff],
	    a = bytes.fromBits(hex.toBits(curve.a.toString(16))),
	    b = bytes.fromBits(hex.toBits(curve.b.toString(16))),
	    xg = bytes.fromBits(hex.toBits(curve.G.x.toString(16))),
	    yg = bytes.fromBits(hex.toBits(curve.G.y.toString(16))),
	    za = entla.concat(ida, a, b, xg, yg, x, y);
	//dumpb(za, "ZA_Src");
	//dumpb(bytes.fromBits(sjcl.hash.sm3.hash(bytes.toBits(za))), "ZA");
	return bytes.fromBits(hash(bytes.toBits(za)));
    /* return hash("1234567812345678"); */  
  },
  
  /**
   * @param {bitArray} za
   * @param {int} klen expected bitLength of output
   * @param {sjcl.hash.xxx.hash} hash function
   * @param {int} bitLength of the output of hash function
   */
  KDF : function(za, klen, hash, v)
  {
	v = v || 256;

	var bytes = sjcl.codec.bytes,
	    ct = 1, 
	    blocks = Math.floor((klen+v-1)/v), 
		ha = [];

	if (klen >= 0xffffffff*v)
	{
		throw new Error("message is too long");
	}
	
    for (var i=0; i<blocks; i++)
	{
		ha = ha.concat(hash(sjcl.bitArray.concat(za, [ct++ & 0xffffffff])));
	}
	return sjcl.bitArray.clamp(ha, klen);
  }
};

/** sm2 publicKey.
 * @constructor
 * @augments sjcl.sm2.basicKey.publicKey
 */
sjcl.ecc.sm2.publicKey = function (curve, point) {
  sjcl.ecc.basicKey.publicKey.apply(this, arguments);
};

/** specific functions for sm2 publicKey. */
sjcl.ecc.sm2.publicKey.prototype = {
  
  /** SM2 verify function
   * @param {string | bitArray} msg message to verify.
   * @param {bitArray} rs signature bitArray.
   * @param {string | bitArray} ida UserID
   * @param {sjcl.hash.xxx.hash} hash hash function 
   */
  verify: function(msg, rs, ida, hash) {
	ida = ida || "1234567812345678";
	hash = hash || sjcl.hash.sm3.hash;
	var bytes = sjcl.codec.bytes,
	string = sjcl.codec.utf8String;
	/* Check whether r and s are in range [1, n-1] */
	var w = sjcl.bitArray,
        n = this._curve.r,
        l = this._curveBitLength,
        r = sjcl.bn.fromBits(w.bitSlice(rs,0,l)),
        s = sjcl.bn.fromBits(w.bitSlice(rs,l,2*l));
     if (r.equals(0) || r.greaterEquals(n) || s.equals(0) || s.greaterEquals(n)) {
      return false;
    }
	if (typeof ida === 'string')
	  ida = bytes.fromBits(string.toBits(ida));
	
    var t = r.add(s).mod(n);
	if (t.equals(0)) {
		return false;
	}
	
	//var x = this._curve.G.mult(s).add(this._point.mult(t)).x;//wr???
	var x = this._curve.G.mult2(s, t, this._point).x;
		
    if (typeof msg === 'string')
	  msg  = bytes.fromBits(string.toBits(msg));
    var za = sjcl.ecc.sm2.HM(this._curve, this._point.x, this._point.y, ida, hash);
	if (sjcl.bitArray.bitLength(za) > this._curveBitLength) {
      za = sjcl.bitArray.clamp(za, this._curveBitLength);
    }
	var e  = sjcl.bn.fromBits(hash(bytes.toBits(za.concat(msg))));

    return e.add(x).mod(n).equals(r);
  },
  /**
   * @param {string | bitArray} msg message to encrypt
   * @param {sjcl.hash.xxx.hash} hash function to be used
   * @param {int} v bitLength of the output of hash function
   */
  encrypt: function(msg, hash, v) {
	hash = hash || sjcl.hash.sm3.hash;
	v = v || 256;
	
	var hex = sjcl.codec.hex,
	    bytes = sjcl.codec.bytes,
		string = sjcl.codec.utf8String;
	
	if (typeof msg === "string") {
		msg = string.toBits(msg)
	}
	var k, c1, c2 = [], za, t, i, nz = 0;

    do {
      k = sjcl.bn.random(this._curve.r.sub(1)).add(1);
	  //k = new sjcl.bn(1);
      c1 = this._curve.G.mult(k);
	  za = this._point.mult(k);/* [x2,y2] */
	  t = sjcl.ecc.sm2.KDF(za.toBits(), sjcl.bitArray.bitLength(msg), hash, v);
      for (i=0; i<t.length; i++)
		nz |= t[i];
	} while (!nz);
	
	for (i=0; i<msg.length; i++) {
		c2[i] = msg[i] ^ t[i];
	}
	za = sjcl.bitArray.concat(za.x.toBits(), sjcl.bitArray.concat(msg, za.y.toBits()));
	var c3 = hash(za);
	return sjcl.bitArray.concat(c1.toBits(), sjcl.bitArray.concat(c3, c2));
  },
  getType: function() {
    return "sm2";
  }
};

/** ecdsa secretKey
* @constructor
* @augments sjcl.ecc.basicKey.publicKey
*/
sjcl.ecc.sm2.secretKey = function (curve, exponent) {
  sjcl.ecc.basicKey.secretKey.apply(this, arguments);
  /* weir007: precompute (1+d)^-1 */
  //Object.assign(this, {'_ivda1':this._exponent.add(1).inverseMod(this._curve.r)});
  this._ivda1 = this._exponent.add(1).inverseMod(this._curve.r);
};

/** specific functions for ecdsa secretKey. */
sjcl.ecc.sm2.secretKey.prototype = {
  /** SM2 sign function
   * @param {string | bitArray} msg message to sign.
   * @param {bitArray} rs signature bitArray.
   * @param {string | bitArray} ida UserID
   * @param {sjcl.hash.xxx.hash} hash hash function 
   */
  sign: function(msg, ida, hash) {
	ida = ida || "1234567812345678";
	hash = hash || sjcl.hash.sm3.hash;
	
	if (typeof msg === 'string')
	  msg  = sjcl.codec.bytes.fromBits(sjcl.codec.utf8String.toBits(msg));
	//if (!this._pk) this._pk = this._curve.G.mult(this._exponent);
	if (!this._pk) this._pk = this._curve.G.multG(this._exponent);
	var za = sjcl.ecc.sm2.HM(this._curve, this._pk.x, this._pk.y, ida, hash);
    if (sjcl.bitArray.bitLength(za) > this._curveBitLength) {
      za = sjcl.bitArray.clamp(za, this._curveBitLength);
    }
    var e  = hash(sjcl.codec.bytes.toBits(za.concat(msg))),
	    n  = this._curve.r, k, x, r, s, k, l = n.bitLength();
		
	do {
	  k  = sjcl.bn.random(n.sub(1)).add(1);
	  //x  = this._curve.G.mult(k).x.mod(n);
	  x  = this._curve.G.multG(k).x.mod(n);
	  r  = sjcl.bn.fromBits(e).add(x).mod(n);
	  k  = k.add(r).mod(n);
	} while(r.equals(0) || k.equals(0));
	//s  = this._exponent.add(1).inverseMod(n).mul(k.mod(n)).sub(r).mod(n);
	s  = this._ivda1.mul(k).sub(r).mod(n);
	if (s.equals(0)) return null;
    
	return sjcl.bitArray.concat(r.toBits(l), s.toBits(l));
  },
  /** SM2 decryption
   * @param {bitArray} cipher ciphertext to decrypt
   * @param {sjcl.hash.xxx.hash} hash hash function 
   * @param {int} v bitLength of the output of hash function
   */
  decrypt: function(cipher, hash, v) {
	hash = hash || sjcl.hash.sm3.hash;
	v = v || 256;

	var klen = sjcl.bitArray.bitLength(cipher) - 768,
	    c1 = sjcl.bitArray.clamp(cipher, 512),
		c3 = sjcl.bitArray.bitSlice(cipher, 512, 768),
		c2 = sjcl.bitArray.bitSlice(cipher, 768),
		i;
		
	c1 = this._curve.fromBits(c1);
	var za = c1.mult(this._exponent),
	    t  = sjcl.ecc.sm2.KDF(za.toBits(), klen, hash, v),
	    nz = 0,
		m = [],
		h;
	
	for (i=0; i<t.length; i++) {
		nz |= t[i];
		m[i] = c2[i] ^ t[i]
	}
	if (!nz) {
	  throw new Error("KDF: encryption key is zero!");
	}
	
	h = hash(sjcl.bitArray.concat(za.x.toBits(), sjcl.bitArray.concat(m, za.y.toBits())));
	
	nz = 0;
	for (i=0; i<h.length; i++)
	  nz |= h[i] ^ c3[i];
    if (nz) {
	  throw new Error("Decryption failed!");
	}
	return m;
  },

  getType: function() {
    return "sm2";
  }
};