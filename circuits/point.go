package circuits

import (
	"fmt"
	"math/big"

	"github.com/consensys/gnark/frontend"
	"github.com/consensys/gnark/std/algebra/algopts"
	"github.com/consensys/gnark/std/math/emulated"
)

// New returns a new [Curve] instance over the base field Base and scalar field
// Scalars defined by the curve parameters params. It returns an error if
// initialising the field emulation fails (for example, when the native field is
// too small) or when the curve parameters are incompatible with the fields.
func New[Base, Scalars emulated.FieldParams](api frontend.API, params CurveParams) (*Curve[Base, Scalars], error) {
	ba, err := emulated.NewField[Base](api)
	if err != nil {
		return nil, fmt.Errorf("new base api: %w", err)
	}
	sa, err := emulated.NewField[Scalars](api)
	if err != nil {
		return nil, fmt.Errorf("new scalar api: %w", err)
	}
	emuGm := make([]AffinePoint[Base], len(params.Gm))
	for i, v := range params.Gm {
		emuGm[i] = AffinePoint[Base]{emulated.ValueOf[Base](v[0]), emulated.ValueOf[Base](v[1])}
	}
	Gx := emulated.ValueOf[Base](params.Gx)
	Gy := emulated.ValueOf[Base](params.Gy)
	var eigenvalue *emulated.Element[Scalars]
	var thirdRootOne *emulated.Element[Base]
	if params.Eigenvalue != nil && params.ThirdRootOne != nil {
		eigenvalue = sa.NewElement(params.Eigenvalue)
		thirdRootOne = ba.NewElement(params.ThirdRootOne)
	}
	return &Curve[Base, Scalars]{
		params:    params,
		api:       api,
		baseApi:   ba,
		scalarApi: sa,
		g: AffinePoint[Base]{
			X: Gx,
			Y: Gy,
		},
		gm:           emuGm,
		a:            emulated.ValueOf[Base](params.A),
		b:            emulated.ValueOf[Base](params.B),
		addA:         params.A.Cmp(big.NewInt(0)) != 0,
		eigenvalue:   eigenvalue,
		thirdRootOne: thirdRootOne,
	}, nil
}

// Curve is an initialised curve which allows performing group operations.
type Curve[Base, Scalars emulated.FieldParams] struct {
	// params is the parameters of the curve
	params CurveParams
	// api is the native api, we construct it ourselves to be sure
	api frontend.API
	// baseApi is the api for point operations
	baseApi *emulated.Field[Base]
	// scalarApi is the api for scalar operations
	scalarApi *emulated.Field[Scalars]

	// g is the generator (base point) of the curve.
	g AffinePoint[Base]

	// gm are the pre-computed doubles the generator (base point) of the curve.
	gm []AffinePoint[Base]

	a            emulated.Element[Base]
	b            emulated.Element[Base]
	addA         bool
	eigenvalue   *emulated.Element[Scalars]
	thirdRootOne *emulated.Element[Base]
}

// Generator returns the base point of the curve. The method does not copy and
// modifying the returned element leads to undefined behaviour!
func (c *Curve[B, S]) Generator() *AffinePoint[B] {
	return &c.g
}

// GeneratorMultiples returns the pre-computed doubles of the base point of the curve. The method does not copy and
// modifying the returned element leads to undefined behaviour!
func (c *Curve[B, S]) GeneratorMultiples() []AffinePoint[B] {
	return c.gm
}

// AffinePoint represents a point on the elliptic curve. We do not check that
// the point is actually on the curve.
//
// Point (0,0) represents point at the infinity. This representation is
// compatible with the EVM representations of points at infinity.
type AffinePoint[Base emulated.FieldParams] struct {
	X, Y emulated.Element[Base]
}

// Neg returns an inverse of p. It doesn't modify p.
func (c *Curve[B, S]) Neg(p *AffinePoint[B]) *AffinePoint[B] {
	return &AffinePoint[B]{
		X: p.X,
		Y: *c.baseApi.Neg(&p.Y),
	}
}

// AssertIsEqual asserts that p and q are the same point.
func (c *Curve[B, S]) AssertIsEqual(p, q *AffinePoint[B]) {
	c.baseApi.AssertIsEqual(&p.X, &q.X)
	c.baseApi.AssertIsEqual(&p.Y, &q.Y)
}

// add adds p and q and returns it. It doesn't modify p nor q.
//
// ⚠️  p must be different than q and -q, and both nonzero.
//
// It uses incomplete formulas in affine coordinates.
func (c *Curve[B, S]) add(p, q *AffinePoint[B]) *AffinePoint[B] {
	mone := c.baseApi.NewElement(-1)
	// compute λ = (q.y-p.y)/(q.x-p.x)
	qypy := c.baseApi.Sub(&q.Y, &p.Y)
	qxpx := c.baseApi.Sub(&q.X, &p.X)
	λ := c.baseApi.Div(qypy, qxpx)

	// xr = λ²-p.x-q.x
	xr := c.baseApi.Eval([][]*emulated.Element[B]{{λ, λ}, {mone, c.baseApi.Add(&p.X, &q.X)}}, []int{1, 1})

	// p.y = λ(p.x-r.x) - p.y
	yr := c.baseApi.Eval([][]*emulated.Element[B]{{λ, c.baseApi.Sub(&p.X, xr)}, {mone, &p.Y}}, []int{1, 1})

	return &AffinePoint[B]{
		X: *c.baseApi.Reduce(xr),
		Y: *c.baseApi.Reduce(yr),
	}
}

// AddUnified adds p and q and returns it. It doesn't modify p nor q.
//
// ✅ p can be equal to q, and either or both can be (0,0).
// (0,0) is not on the curve but we conventionally take it as the
// neutral/infinity point as per the [EVM].
//
// It uses the unified formulas of Brier and Joye ([[BriJoy02]] (Corollary 1)).
//
// [BriJoy02]: https://link.springer.com/content/pdf/10.1007/3-540-45664-3_24.pdf
// [EVM]: https://ethereum.github.io/yellowpaper/paper.pdf
func (c *Curve[B, S]) AddUnified(p, q *AffinePoint[B]) *AffinePoint[B] {

	// selector1 = 1 when p is (0,0) and 0 otherwise
	selector1 := c.api.And(c.baseApi.IsZero(&p.X), c.baseApi.IsZero(&p.Y))
	// selector2 = 1 when q is (0,0) and 0 otherwise
	selector2 := c.api.And(c.baseApi.IsZero(&q.X), c.baseApi.IsZero(&q.Y))

	// λ = ((p.x+q.x)² - p.x*q.x + a)/(p.y + q.y)
	pxqx := c.baseApi.MulMod(&p.X, &q.X)
	pxplusqx := c.baseApi.Add(&p.X, &q.X)
	num := c.baseApi.MulMod(pxplusqx, pxplusqx)
	num = c.baseApi.Sub(num, pxqx)
	if c.addA {
		num = c.baseApi.Add(num, &c.a)
	}
	denum := c.baseApi.Add(&p.Y, &q.Y)
	// if p.y + q.y = 0, assign dummy 1 to denum and continue
	selector3 := c.baseApi.IsZero(denum)
	denum = c.baseApi.Select(selector3, c.baseApi.One(), denum)
	λ := c.baseApi.Div(num, denum)

	// x = λ^2 - p.x - q.x
	xr := c.baseApi.MulMod(λ, λ)
	xr = c.baseApi.Sub(xr, pxplusqx)

	// y = λ(p.x - xr) - p.y
	yr := c.baseApi.Sub(&p.X, xr)
	yr = c.baseApi.MulMod(yr, λ)
	yr = c.baseApi.Sub(yr, &p.Y)
	result := AffinePoint[B]{
		X: *c.baseApi.Reduce(xr),
		Y: *c.baseApi.Reduce(yr),
	}

	zero := c.baseApi.Zero()
	infinity := AffinePoint[B]{X: *zero, Y: *zero}
	// if p=(0,0) return q
	result = *c.Select(selector1, q, &result)
	// if q=(0,0) return p
	result = *c.Select(selector2, p, &result)
	// if p.y + q.y = 0, return (0, 0)
	result = *c.Select(selector3, &infinity, &result)

	return &result
}

// double doubles p and return it. It doesn't modify p.
//
// ⚠️  p.Y must be nonzero.
//
// It uses affine coordinates.
func (c *Curve[B, S]) double(p *AffinePoint[B]) *AffinePoint[B] {

	mone := c.baseApi.NewElement(-1)
	// compute λ = (3p.x²+a)/2*p.y, here we assume a=0 (j invariant 0 curve)
	xx3a := c.baseApi.MulMod(&p.X, &p.X)
	xx3a = c.baseApi.MulConst(xx3a, big.NewInt(3))
	if c.addA {
		xx3a = c.baseApi.Add(xx3a, &c.a)
	}
	y2 := c.baseApi.MulConst(&p.Y, big.NewInt(2))
	λ := c.baseApi.Div(xx3a, y2)

	// xr = λ²-2p.x
	xr := c.baseApi.Eval([][]*emulated.Element[B]{{λ, λ}, {mone, &p.X}}, []int{1, 2})

	// yr = λ(p-xr) - p.y
	yr := c.baseApi.Eval([][]*emulated.Element[B]{{λ, c.baseApi.Sub(&p.X, xr)}, {mone, &p.Y}}, []int{1, 1})

	return &AffinePoint[B]{
		X: *c.baseApi.Reduce(xr),
		Y: *c.baseApi.Reduce(yr),
	}
}

// triple triples p and return it. It follows [ELM03] (Section 3.1).
// Saves the computation of the y coordinate of 2p as it is used only in the computation of λ2,
// which can be computed as
//
//	λ2 = -λ1-2*p.y/(x2-p.x)
//
// instead. It doesn't modify p.
//
// ⚠️  p.Y must be nonzero.
//
// [ELM03]: https://arxiv.org/pdf/math/0208038.pdf
func (c *Curve[B, S]) triple(p *AffinePoint[B]) *AffinePoint[B] {

	mone := c.baseApi.NewElement(-1)
	// compute λ1 = (3p.x²+a)/2p.y, here we assume a=0 (j invariant 0 curve)
	xx := c.baseApi.MulMod(&p.X, &p.X)
	xx = c.baseApi.MulConst(xx, big.NewInt(3))
	if c.addA {
		xx = c.baseApi.Add(xx, &c.a)
	}
	y2 := c.baseApi.MulConst(&p.Y, big.NewInt(2))
	λ1 := c.baseApi.Div(xx, y2)

	// xr = λ1²-2p.x
	x2 := c.baseApi.Eval([][]*emulated.Element[B]{{λ1, λ1}, {mone, &p.X}}, []int{1, 2})

	// omit y2 computation, and
	// compute λ2 = 2p.y/(x2 − p.x) − λ1.
	x1x2 := c.baseApi.Sub(&p.X, x2)
	λ2 := c.baseApi.Div(y2, x1x2)
	λ2 = c.baseApi.Sub(λ2, λ1)

	// xr = λ²-p.x-x2
	xr := c.baseApi.Eval([][]*emulated.Element[B]{{λ2, λ2}, {mone, &p.X}, {mone, x2}}, []int{1, 1, 1})

	// yr = λ(p.x-xr) - p.y
	yr := c.baseApi.Eval([][]*emulated.Element[B]{{λ2, c.baseApi.Sub(&p.X, xr)}, {mone, &p.Y}}, []int{1, 1})

	return &AffinePoint[B]{
		X: *c.baseApi.Reduce(xr),
		Y: *c.baseApi.Reduce(yr),
	}
}

// doubleAndAdd computes 2p+q as (p+q)+p. It follows [ELM03] (Section 3.1)
// Saves the computation of the y coordinate of p+q as it is used only in the computation of λ2,
// which can be computed as
//
//	λ2 = -λ1-2*p.y/(x2-p.x)
//
// instead. It doesn't modify p nor q.
//
// ⚠️  p must be different than q and -q, and both nonzero.
//
// [ELM03]: https://arxiv.org/pdf/math/0208038.pdf
func (c *Curve[B, S]) doubleAndAdd(p, q *AffinePoint[B]) *AffinePoint[B] {

	mone := c.baseApi.NewElement(-1)
	// compute λ1 = (q.y-p.y)/(q.x-p.x)
	yqyp := c.baseApi.Sub(&q.Y, &p.Y)
	xpn := c.baseApi.Neg(&p.X)
	xqxp := c.baseApi.Add(&q.X, xpn)
	λ1 := c.baseApi.Div(yqyp, xqxp)

	// compute x2 = λ1²-p.x-q.x
	x2 := c.baseApi.Eval([][]*emulated.Element[B]{{λ1, λ1}, {mone, c.baseApi.Add(&p.X, &q.X)}}, []int{1, 1})

	// omit y2 computation

	// compute -λ2 = λ1+2*p.y/(x2-p.x)
	ypyp := c.baseApi.MulConst(&p.Y, big.NewInt(2))
	x2xp := c.baseApi.Add(x2, xpn)
	λ2 := c.baseApi.Div(ypyp, x2xp)
	λ2 = c.baseApi.Add(λ1, λ2)

	// compute x3 = (-λ2)²-p.x-x2
	x3 := c.baseApi.Eval([][]*emulated.Element[B]{{λ2, λ2}, {mone, &p.X}, {mone, x2}}, []int{1, 1, 1})

	// compute y3 = -λ2*(x3 - p.x)-p.y
	y3 := c.baseApi.Eval([][]*emulated.Element[B]{{λ2, c.baseApi.Add(x3, xpn)}, {mone, &p.Y}}, []int{1, 1})

	return &AffinePoint[B]{
		X: *c.baseApi.Reduce(x3),
		Y: *c.baseApi.Reduce(y3),
	}

}

// scalarMulFakeGLV computes [s]Q and returns it. It doesn't modify Q nor s.
//
// ⚠️  The scalar s must be nonzero and the point Q different from (0,0) unless [algopts.WithCompleteArithmetic] is set.
// (0,0) is not on the curve but we conventionally take it as the
// neutral/infinity point as per the [EVM].
//
// [EVM]: https://ethereum.github.io/yellowpaper/paper.pdf
func (c *Curve[B, S]) scalarMulFakeGLV(Q *AffinePoint[B], s *emulated.Element[S], opts ...algopts.AlgebraOption) *AffinePoint[B] {
	cfg, err := algopts.NewConfig(opts...)
	if err != nil {
		panic(err)
	}

	var selector1 frontend.Variable
	_s := s
	if cfg.CompleteArithmetic {
		selector1 = c.scalarApi.IsZero(s)
		_s = c.scalarApi.Select(selector1, c.scalarApi.One(), s)
	}

	// First we find the sub-salars s1, s2 s.t. s1 + s2*s = 0 mod r and s1, s2 < sqrt(r).
	sd, err := c.scalarApi.NewHint(halfGCD, 2, _s)
	if err != nil {
		panic(fmt.Sprintf("halfGCD hint: %v", err))
	}
	s1, s2 := sd[0], sd[1]
	// s2 can be negative. If so, we return in the halfGCD hint -s2
	// and here compute _s2 = -s2 mod r
	sign, err := c.scalarApi.NewHintWithNativeOutput(halfGCDSigns, 1, _s)
	if err != nil {
		panic(fmt.Sprintf("halfGCDSigns hint: %v", err))
	}
	_s2 := c.scalarApi.Select(sign[0], c.scalarApi.Neg(s2), s2)
	// We check that s1 + s*_s2 == 0 mod r
	c.scalarApi.AssertIsEqual(
		c.scalarApi.Add(s1, c.scalarApi.Mul(_s, _s2)),
		c.scalarApi.Zero(),
	)
	// A malicious hint can provide s1=s2=0 mod r
	// So we check that _s2 is non-zero otherwise [0]([s]Q = ∀R) is always true
	c.api.AssertIsEqual(c.scalarApi.IsZero(_s2), 0)

	// Then we compute the hinted scalar mul R = [s]Q
	// Q coordinates are in Fp and the scalar s in Fr
	// we decompose Q.X, Q.Y, s into limbs and recompose them in the hint.
	var inps []frontend.Variable
	inps = append(inps, Q.X.Limbs...)
	inps = append(inps, Q.Y.Limbs...)
	inps = append(inps, s.Limbs...)
	R, err := c.baseApi.NewHintWithNativeInput(scalarMulHint, 2, inps...)
	if err != nil {
		panic(fmt.Sprintf("scalar mul hint: %v", err))
	}
	r0, r1 := R[0], R[1]

	var selector2 frontend.Variable
	one := c.baseApi.One()
	dummy := &AffinePoint[B]{X: *one, Y: *one}
	addFn := c.add
	if cfg.CompleteArithmetic {
		addFn = c.AddUnified
		// if Q=(0,0) we assign a dummy (1,1) to Q and R and continue
		selector2 = c.api.And(c.baseApi.IsZero(&Q.X), c.baseApi.IsZero(&Q.Y))
		Q = c.Select(selector2, dummy, Q)
		r0 = c.baseApi.Select(selector2, c.baseApi.Zero(), r0)
		r1 = c.baseApi.Select(selector2, &dummy.Y, r1)
	}

	var st S
	nbits := (st.Modulus().BitLen() + 1) / 2
	s1bits := c.scalarApi.ToBits(s1)
	s2bits := c.scalarApi.ToBits(s2)

	// Precomputations:
	// 		tableQ[0] = -Q
	//   	tableQ[1] = Q
	// 		tableQ[2] = [3]Q
	// 		tableR[0] = -R or R if s2 is negative
	//   	tableR[1] = R or -R if s2 is negative
	// 		tableR[2] = [3]R or [-3]R if s2 is negative
	var tableQ, tableR [3]*AffinePoint[B]
	tableQ[1] = Q
	tableQ[0] = c.Neg(Q)
	tableQ[2] = c.triple(tableQ[1])
	tableR[1] = &AffinePoint[B]{
		X: *r0,
		Y: *c.baseApi.Select(sign[0], c.baseApi.Neg(r1), r1),
	}
	tableR[0] = c.Neg(tableR[1])
	if cfg.CompleteArithmetic {
		tableR[2] = c.AddUnified(tableR[1], tableR[1])
		tableR[2] = c.AddUnified(tableR[2], tableR[1])
	} else {
		tableR[2] = c.triple(tableR[1])
	}

	// We should start the accumulator by the infinity point, but since affine
	// formulae are incomplete we suppose that the first bits of the
	// sub-scalars s1 and s2 are 1, and set:
	// 		Acc = Q + R
	Acc := c.add(tableQ[1], tableR[1])

	// At each iteration we need to compute:
	// 		[2]Acc ± Q ± R.
	// We can compute [2]Acc and look up the (precomputed) point P from:
	// 		B1 = Q+R
	// 		B2 = -Q-R
	// 		B3 = Q-R
	// 		B4 = -Q+R
	//
	// If we extend this by merging two iterations, we need to look up P and P'
	// both from {B1, B2, B3, B4} and compute:
	// 		[2]([2]Acc+P)+P' = [4]Acc + T
	// where T = [2]P+P'. So at each (merged) iteration, we can compute [4]Acc
	// and look up T from the precomputed list of points:
	//
	// T = [3](Q + R)
	// P = B1 and P' = B1
	T1 := c.add(tableQ[2], tableR[2])
	// T = Q + R
	// P = B1 and P' = B2
	T2 := Acc
	// T = [3]Q + R
	// P = B1 and P' = B3
	T3 := c.add(tableQ[2], tableR[1])
	// T = Q + [3]R
	// P = B1 and P' = B4
	T4 := c.add(tableQ[1], tableR[2])
	// T  = -Q - R
	// P = B2 and P' = B1
	T5 := c.Neg(T2)
	// T  = -[3](Q + R)
	// P = B2 and P' = B2
	T6 := c.Neg(T1)
	// T = -Q - [3]R
	// P = B2 and P' = B3
	T7 := c.Neg(T4)
	// T = -[3]Q - R
	// P = B2 and P' = B4
	T8 := c.Neg(T3)
	// T = [3]Q - R
	// P = B3 and P' = B1
	T9 := c.add(tableQ[2], tableR[0])
	// T = Q - [3]R
	// P = B3 and P' = B2
	T11 := c.Neg(tableR[2])
	T10 := c.add(tableQ[1], T11)
	// T = [3](Q - R)
	// P = B3 and P' = B3
	T11 = c.add(tableQ[2], T11)
	// T = -R + Q
	// P = B3 and P' = B4
	T12 := c.add(tableR[0], tableQ[1])
	// T = [3]R - Q
	// P = B4 and P' = B1
	T13 := c.Neg(T10)
	// T = R - [3]Q
	// P = B4 and P' = B2
	T14 := c.Neg(T9)
	// T = R - Q
	// P = B4 and P' = B3
	T15 := c.Neg(T12)
	// T = [3](R - Q)
	// P = B4 and P' = B4
	T16 := c.Neg(T11)
	// note that half of these points are negatives of the other half,
	// hence have the same X coordinates.

	// When nbits is even, we need to handle the first iteration separately
	if nbits%2 == 0 {
		// Acc = [2]Acc ± Q ± R
		T := &AffinePoint[B]{
			X: *c.baseApi.Select(c.api.Xor(s1bits[nbits-1], s2bits[nbits-1]), &T12.X, &T5.X),
			Y: *c.baseApi.Lookup2(s1bits[nbits-1], s2bits[nbits-1], &T5.Y, &T12.Y, &T15.Y, &T2.Y),
		}
		// We don't use doubleAndAdd here as it would involve edge cases
		// when bits are 00 (T==-Acc) or 11 (T==Acc).
		Acc = c.double(Acc)
		Acc = c.add(Acc, T)
	} else {
		// when nbits is odd we start the main loop at normally nbits - 1
		nbits++
	}
	for i := nbits - 2; i > 2; i -= 2 {
		// selectorY takes values in [0,15]
		selectorY := c.api.Add(
			s1bits[i],
			c.api.Mul(s2bits[i], 2),
			c.api.Mul(s1bits[i-1], 4),
			c.api.Mul(s2bits[i-1], 8),
		)
		// selectorX takes values in [0,7] s.t.:
		// 		- when selectorY < 8: selectorX = selectorY
		// 		- when selectorY >= 8: selectorX = 15 - selectorY
		selectorX := c.api.Add(
			c.api.Mul(selectorY, c.api.Sub(1, c.api.Mul(s2bits[i-1], 2))),
			c.api.Mul(s2bits[i-1], 15),
		)
		// Bi.Y are distincts so we need a 16-to-1 multiplexer,
		// but only half of the Bi.X are distinct so we need a 8-to-1.
		T := &AffinePoint[B]{
			X: *c.baseApi.Mux(selectorX,
				&T6.X, &T10.X, &T14.X, &T2.X, &T7.X, &T11.X, &T15.X, &T3.X,
			),
			Y: *c.baseApi.Mux(selectorY,
				&T6.Y, &T10.Y, &T14.Y, &T2.Y, &T7.Y, &T11.Y, &T15.Y, &T3.Y,
				&T8.Y, &T12.Y, &T16.Y, &T4.Y, &T5.Y, &T9.Y, &T13.Y, &T1.Y,
			),
		}
		// Acc = [4]Acc + T
		Acc = c.double(Acc)
		Acc = c.doubleAndAdd(Acc, T)
	}

	// i = 2
	// we isolate the last iteration to avoid falling into incomplete additions
	//
	// selectorY takes values in [0,15]
	selectorY := c.api.Add(
		s1bits[2],
		c.api.Mul(s2bits[2], 2),
		c.api.Mul(s1bits[1], 4),
		c.api.Mul(s2bits[1], 8),
	)
	// selectorX takes values in [0,7] s.t.:
	// 		- when selectorY < 8: selectorX = selectorY
	// 		- when selectorY >= 8: selectorX = 15 - selectorY
	selectorX := c.api.Add(
		c.api.Mul(selectorY, c.api.Sub(1, c.api.Mul(s2bits[1], 2))),
		c.api.Mul(s2bits[1], 15),
	)
	// Bi.Y are distincts so we need a 16-to-1 multiplexer,
	// but only half of the Bi.X are distinct so we need a 8-to-1.
	T := &AffinePoint[B]{
		X: *c.baseApi.Mux(selectorX,
			&T6.X, &T10.X, &T14.X, &T2.X, &T7.X, &T11.X, &T15.X, &T3.X,
		),
		Y: *c.baseApi.Mux(selectorY,
			&T6.Y, &T10.Y, &T14.Y, &T2.Y, &T7.Y, &T11.Y, &T15.Y, &T3.Y,
			&T8.Y, &T12.Y, &T16.Y, &T4.Y, &T5.Y, &T9.Y, &T13.Y, &T1.Y,
		),
	}
	// to avoid incomplete additions we add [3]R to the precomputed T before computing [4]Acc+T
	// 		Acc = [4]Acc + T + [3]R
	T = c.add(T, tableR[2])
	Acc = c.double(Acc)
	Acc = c.doubleAndAdd(Acc, T)

	// i = 0
	// subtract Q and R if the first bits are 0.
	// When cfg.CompleteArithmetic is set, we use AddUnified instead of Add.
	// This means when s=0 then Acc=(0,0) because AddUnified(Q, -Q) = (0,0).
	tableQ[0] = addFn(tableQ[0], Acc)
	Acc = c.Select(s1bits[0], Acc, tableQ[0])
	tableR[0] = addFn(tableR[0], Acc)
	Acc = c.Select(s2bits[0], Acc, tableR[0])

	if cfg.CompleteArithmetic {
		Acc = c.Select(c.api.Or(selector1, selector2), tableR[2], Acc)
	}
	// we added [3]R at the last iteration so the result should be
	// 		Acc = [s1]Q + [s2]R + [3]R
	// 		    = [s1]Q + [s2*s]Q + [3]R
	// 		    = [s1+s2*s]Q + [3]R
	// 		    = [0]Q + [3]R
	// 		    = [3]R
	c.AssertIsEqual(Acc, tableR[2])

	return &AffinePoint[B]{
		X: *R[0],
		Y: *R[1],
	}
}

// Select selects between p and q given the selector b. If b == 1, then returns
// p and q otherwise.
func (c *Curve[B, S]) Select(b frontend.Variable, p, q *AffinePoint[B]) *AffinePoint[B] {
	x := c.baseApi.Select(b, &p.X, &q.X)
	y := c.baseApi.Select(b, &p.Y, &q.Y)
	return &AffinePoint[B]{
		X: *x,
		Y: *y,
	}
}

// Lookup2 performs a 2-bit lookup between i0, i1, i2, i3 based on bits b0
// and b1. Returns:
//   - i0 if b0=0 and b1=0,
//   - i1 if b0=1 and b1=0,
//   - i2 if b0=0 and b1=1,
//   - i3 if b0=1 and b1=1.
func (c *Curve[B, S]) Lookup2(b0, b1 frontend.Variable, i0, i1, i2, i3 *AffinePoint[B]) *AffinePoint[B] {
	x := c.baseApi.Lookup2(b0, b1, &i0.X, &i1.X, &i2.X, &i3.X)
	y := c.baseApi.Lookup2(b0, b1, &i0.Y, &i1.Y, &i2.Y, &i3.Y)
	return &AffinePoint[B]{
		X: *x,
		Y: *y,
	}
}

// Mux performs a lookup from the inputs and returns inputs[sel]. It is most
// efficient for power of two lengths of the inputs, but works for any number of
// inputs.
func (c *Curve[B, S]) Mux(sel frontend.Variable, inputs ...*AffinePoint[B]) *AffinePoint[B] {
	xs := make([]*emulated.Element[B], len(inputs))
	ys := make([]*emulated.Element[B], len(inputs))
	for i := range inputs {
		xs[i] = &inputs[i].X
		ys[i] = &inputs[i].Y
	}
	return &AffinePoint[B]{
		X: *c.baseApi.Mux(sel, xs...),
		Y: *c.baseApi.Mux(sel, ys...),
	}
}

// scalarMulJoye computes [s]p and returns it. It doesn't modify p nor s.
// This function doesn't check that the p is on the curve. See AssertIsOnCurve.
//
// ⚠️  p must not be (0,0) and s must not be 0, unless [algopts.WithCompleteArithmetic] option is set.
// (0,0) is not on the curve but we conventionally take it as the
// neutral/infinity point as per the [EVM].
//
// It computes the right-to-left variable-base double-and-add algorithm ([Joye07], Alg.1).
//
// Since we use incomplete formulas for the addition law, we need to start with
// a non-zero accumulator point (R0). To do this, we skip the LSB (bit at
// position 0) and proceed assuming it was 1. At the end, we conditionally
// subtract the initial value (p) if LSB is 1. We also handle the bits at
// positions 1 and n-1 outside of the loop to optimize the number of
// constraints using [ELM03] (Section 3.1)
//
// [ELM03]: https://arxiv.org/pdf/math/0208038.pdf
// [EVM]: https://ethereum.github.io/yellowpaper/paper.pdf
// [Joye07]: https://www.iacr.org/archive/ches2007/47270135/47270135.pdf
func (c *Curve[B, S]) scalarMulJoye(p *AffinePoint[B], s *emulated.Element[S], opts ...algopts.AlgebraOption) *AffinePoint[B] {
	cfg, err := algopts.NewConfig(opts...)
	if err != nil {
		panic(fmt.Sprintf("parse opts: %v", err))
	}
	var selector frontend.Variable
	if cfg.CompleteArithmetic {
		// if p=(0,0) we assign a dummy (0,1) to p and continue
		selector = c.api.And(c.baseApi.IsZero(&p.X), c.baseApi.IsZero(&p.Y))
		one := c.baseApi.One()
		p = c.Select(selector, &AffinePoint[B]{X: *one, Y: *one}, p)
	}

	var st S
	sr := c.scalarApi.Reduce(s)
	sBits := c.scalarApi.ToBits(sr)
	n := st.Modulus().BitLen()
	if cfg.NbScalarBits > 2 && cfg.NbScalarBits < n {
		n = cfg.NbScalarBits
	}

	// i = 1
	Rb := c.triple(p)
	R0 := c.Select(sBits[1], Rb, p)
	R1 := c.Select(sBits[1], p, Rb)

	for i := 2; i < n-1; i++ {
		Rb = c.doubleAndAddSelect(sBits[i], R0, R1)
		R0 = c.Select(sBits[i], Rb, R0)
		R1 = c.Select(sBits[i], R1, Rb)
	}

	// i = n-1
	Rb = c.doubleAndAddSelect(sBits[n-1], R0, R1)
	R0 = c.Select(sBits[n-1], Rb, R0)

	// i = 0
	// we use AddUnified instead of Add. This is because:
	// 		- when s=0 then R0=P and AddUnified(P, -P) = (0,0). We return (0,0).
	// 		- when s=1 then R0=P AddUnified(Q, -Q) is well defined. We return R0=P.
	R0 = c.Select(sBits[0], R0, c.AddUnified(R0, c.Neg(p)))

	if cfg.CompleteArithmetic {
		// if p=(0,0), return (0,0)
		zero := c.baseApi.Zero()
		R0 = c.Select(selector, &AffinePoint[B]{X: *zero, Y: *zero}, R0)
	}

	return R0
}

// doubleAndAddSelect is the same as doubleAndAdd but computes either:
//
//	2p+q if b=1 or
//	2q+p if b=0
//
// It first computes the x-coordinate of p+q via the slope(p,q)
// and then based on a Select adds either p or q.
func (c *Curve[B, S]) doubleAndAddSelect(b frontend.Variable, p, q *AffinePoint[B]) *AffinePoint[B] {

	mone := c.baseApi.NewElement(-1)
	// compute λ1 = (q.y-p.y)/(q.x-p.x)
	yqyp := c.baseApi.Sub(&q.Y, &p.Y)
	xqxp := c.baseApi.Sub(&q.X, &p.X)
	λ1 := c.baseApi.Div(yqyp, xqxp)

	// compute x2 = λ1²-p.x-q.x
	x2 := c.baseApi.Eval([][]*emulated.Element[B]{{λ1, λ1}, {mone, &p.X}, {mone, &q.X}}, []int{1, 1, 1})

	// omit y2 computation

	// conditional second addition
	t := c.Select(b, p, q)

	// compute -λ2 = λ1+2*t.y/(x2-t.x)
	ypyp := c.baseApi.MulConst(&t.Y, big.NewInt(2))
	x2xp := c.baseApi.Sub(x2, &t.X)
	λ2 := c.baseApi.Div(ypyp, x2xp)
	λ2 = c.baseApi.Add(λ1, λ2)

	// compute x3 = (-λ2)²-t.x-x2
	x3 := c.baseApi.Eval([][]*emulated.Element[B]{{λ2, λ2}, {mone, &t.X}, {mone, x2}}, []int{1, 1, 1})

	// compute y3 = -λ2*(x3 - t.x)-t.y
	y3 := c.baseApi.Eval([][]*emulated.Element[B]{{λ2, x3}, {mone, λ2, &t.X}, {mone, &t.Y}}, []int{1, 1, 1})

	return &AffinePoint[B]{
		X: *c.baseApi.Reduce(x3),
		Y: *c.baseApi.Reduce(y3),
	}

}
