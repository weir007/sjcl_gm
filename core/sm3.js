/** @fileOverview Javascript SM3 implementation.
 *
 *
 * @author weir007
 */

/**
 * Context for a SM3 operation in progress.
 * @constructor
 */
sjcl.hash.sm3 = function (hash) {
  if (hash) {/* copy the context/state */
    this._h = hash._h.slice(0);/* h is V */
    this._buffer = hash._buffer.slice(0);
    this._length = hash._length;
  } else {
    this.reset();
  }
};

/**
 * Hash a string or an array of words.
 * @static
 * @param {bitArray|String} data the data to hash.
 * @return {bitArray} The hash value, an array of 16 big-endian words.
 */
sjcl.hash.sm3.hash = function (data) {
  return (new sjcl.hash.sm3()).update(data).finalize();
};

sjcl.hash.sm3.prototype = {
  /**
   * The hash's block size, in bits.
   * @constant
   */
  blockSize: 512,
   
  /**
   * Reset the hash state.
   * @return this
   */
  reset:function () {
    this._h = this._init.slice(0);
    this._buffer = [];
    this._length = 0;
    return this;
  },
  
  /**
   * Input several words to the hash.
   * @param {bitArray|String} data the data to hash.
   * @return this
   */
  update: function (data) {
    if (typeof data === "string") {
      data = sjcl.codec.utf8String.toBits(data);
    }
    var i, b = this._buffer = sjcl.bitArray.concat(this._buffer, data),
        ol = this._length,
        nl = this._length = ol + sjcl.bitArray.bitLength(data);
    if (nl > 9007199254740991){
      throw new sjcl.exception.invalid("Cannot hash more than 2^53 - 1 bits");//why
    }

    if (typeof Uint32Array !== 'undefined') {
		var c = new Uint32Array(b);
    	var j = 0;
		for (i = 512+ol - ((512+ol) & 511); i <= nl; i+= 512) {
      	    this._block(c.subarray(16 * j, 16 * (j+1)));
      	    //console.log('V'+j+':'+this._h[0].toString(16)+this._h[1].toString(16)+this._h[2].toString(16)+this._h[3].toString(16)+this._h[4].toString(16)+this._h[5].toString(16)+this._h[6].toString(16)+this._h[7].toString(16));
			j += 1;
    	}
    	b.splice(0, 16 * j);
    } else {
		for (i = 512+ol - ((512+ol) & 511); i <= nl; i+= 512) {
      	    this._block(b.splice(0,16));
      	}
    }
    return this;
  },
  
  /**
   * Complete hashing and output the hash value.
   * @return {bitArray} The hash value, an array of 8 big-endian words.
   */
  finalize:function () {
    var i, b = this._buffer, h = this._h;
    // Round out and push the buthis._FFer
    b = sjcl.bitArray.concat(b, [sjcl.bitArray.partial(1,1)]);
    
    // Round out the buthis._FFer to a multiple of 16 words, less the 2 length words.
    for (i = b.length + 2; i & 15; i++) {
      b.push(0);
    }
    
    // append the length
    b.push(Math.floor(this._length / 0x100000000));
    b.push(this._length | 0);

    while (b.length) {
     this._block(b.splice(0,16));
    }

    this.reset();
    return h;
  },

  /**
   * The SM3 initialization vector.
   * @private
   */
  _init:[0x7380166f,0x4914b2b9,0x172442d7,0xda8a0600,0xa96f30bc,0x163138aa,0xe38dee4d,0xb0fb0e4e],
  
  /**
   * functions to be used in update(_block)
   * _rotl, _P0, _P1, _FF, _GG, _T
   */
  _rotl:function (w, n) {
	return (w<<n) | (w >>> (32-n));
  },
  
  _P0:function (X) {
	return X ^ this._rotl(X, 9) ^ this._rotl(X, 17);
  },

  _P1:function (X) {
	return X ^ this._rotl(X, 15) ^ this._rotl(X, 23);
  },

  _FF:function (X, Y, Z, j) {
	return j >= 0 && j <= 15 ? X ^ Y ^ Z : (X & Y) | (X & Z) | (Y & Z);
  },

  _GG:function (X, Y, Z, j) {
	return j >= 0 && j <= 15 ? X ^ Y ^ Z : (X & Y) | (~X & Z);
  },

  _T:function(j) {
	return j >= 0 && j <= 15 ? 0x79cc4519 : 0x7a879d8a;
  },
  /**
   * Perform one cycle of SM3.
   * @param {Uint32Array|bitArray} w one block of words.
   * @private
   */
  _block:function (w) {  
  var h = this._h,
  W = [],
  M = [];// W'
  for (var i=0; i<16; i++)
	  W[i] = w[i];
  
  for (var j = 16; j < 68; j += 1) {
    W[j] = this._P1(W[j - 16] ^ W[j - 9] ^ this._rotl(W[j - 3], 15)) ^ this._rotl(W[j - 13], 7) ^ W[j - 6];
  }

  // W′[j] = W[j] xor W[j+4]
  for (var j = 0; j < 64; j += 1) {
    M.push(W[j] ^ W[j + 4]);
  }

  var A = h[0], B = h[1], C = h[2], D = h[3], E = h[4], F = h[5], G = h[6], H = h[7];
 
  var SS1, SS2, TT1, TT2;
 
  for (var j = 0; j < 64; j += 1) {
	SS2 = this._rotl(A, 12);
	SS1 = this._rotl((SS2 + E + this._rotl(this._T(j), j & 0x1f)) & 0xffffffff, 7);
	SS2 ^= SS1;
    
    TT1 = this._FF(A, B, C, j) + D + SS2 + M[j];
    TT2 = this._GG(E, F, G, j) + H + SS1 + W[j];

    D = C;
    C = this._rotl(B, 9);
    B = A;
    A = TT1;
    H = G;
    G = this._rotl(F, 19);
    F = E;
    E = this._P0(TT2);
  }

  h[0] ^= A;
  h[1] ^= B;
  h[2] ^= C;
  h[3] ^= D;
  h[4] ^= E;
  h[5] ^= F;
  h[6] ^= G;
  h[7] ^= H;
}
};


