package circuits

import (
	"crypto/elliptic"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc/bn254"
	"github.com/consensys/gnark/std/math/emulated"
)

// CurveParams defines parameters of an elliptic curve in short Weierstrass form
// given by the equation
//
//	Y² = X³ + aX + b
//
// The base point is defined by (Gx, Gy).
type CurveParams struct {
	A            *big.Int      // a in curve equation
	B            *big.Int      // b in curve equation
	Gx           *big.Int      // base point x
	Gy           *big.Int      // base point y
	Gm           [][2]*big.Int // m*base point coords
	Eigenvalue   *big.Int      // endomorphism eigenvalue
	ThirdRootOne *big.Int      // endomorphism image scaler
}

// GetBN254Params returns the curve parameters for the curve BN254 (alt_bn128).
// When initialising new curve, use the base field [emulated.BN254Fp] and scalar
// field [emulated.BN254Fr].
func GetBN254Params() CurveParams {
	_, _, g1aff, _ := bn254.Generators()
	lambda, _ := new(big.Int).SetString("4407920970296243842393367215006156084916469457145843978461", 10)
	omega, _ := new(big.Int).SetString("2203960485148121921418603742825762020974279258880205651966", 10)
	return CurveParams{
		A:            big.NewInt(0),
		B:            big.NewInt(3),
		Gx:           g1aff.X.BigInt(new(big.Int)),
		Gy:           g1aff.Y.BigInt(new(big.Int)),
		Gm:           computeBN254Table(),
		Eigenvalue:   lambda,
		ThirdRootOne: omega,
	}
}

// GetP256Params returns the curve parameters for the curve P-256 (also
// SECP256r1). When initialising new curve, use the base field
// [emulated.P256Fp] and scalar field [emulated.P256Fr].
func GetP256Params() CurveParams {
	pr := elliptic.P256().Params()
	a := new(big.Int).Sub(pr.P, big.NewInt(3))
	return CurveParams{
		A:            a,
		B:            pr.B,
		Gx:           pr.Gx,
		Gy:           pr.Gy,
		Gm:           computeP256Table(),
		Eigenvalue:   nil,
		ThirdRootOne: nil,
	}
}

// GetCurveParams returns suitable curve parameters given the parametric type
// Base as base field. It caches the parameters and modifying the values in the
// parameters struct leads to undefined behaviour.
func GetCurveParams[Base emulated.FieldParams]() CurveParams {
	var t Base
	switch t.Modulus().String() {
	case emulated.BN254Fp{}.Modulus().String():
		return bn254Params
	case emulated.P256Fp{}.Modulus().String():
		return p256Params
	default:
		panic("no stored parameters")
	}
}

var (
	bn254Params CurveParams
	p256Params  CurveParams
)

func init() {
	bn254Params = GetBN254Params()
	p256Params = GetP256Params()
}
