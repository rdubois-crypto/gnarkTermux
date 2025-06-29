package circuits

import (
	"crypto/elliptic"
	"errors"
	"fmt"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark-crypto/ecc/bn254"
	bn_fp "github.com/consensys/gnark-crypto/ecc/bn254/fp"
	"github.com/consensys/gnark/constraint/solver"
	"github.com/consensys/gnark/std/math/emulated"
	"github.com/consensys/gnark/std/math/emulated/emparams"
)

func init() {
	solver.RegisterHint(GetHints()...)
}

func GetHints() []solver.Hint {
	return []solver.Hint{
		decomposeScalarG1Signs,
		decomposeScalarG1Subscalars,
		scalarMulHint,
		halfGCD,
		halfGCDSigns,
	}
}

func decomposeScalarG1Subscalars(mod *big.Int, inputs []*big.Int, outputs []*big.Int) error {
	return emulated.UnwrapHint(inputs, outputs, func(field *big.Int, inputs, outputs []*big.Int) error {
		if len(inputs) != 2 {
			return errors.New("expecting two inputs")
		}
		if len(outputs) != 2 {
			return errors.New("expecting two outputs")
		}
		glvBasis := new(ecc.Lattice)
		ecc.PrecomputeLattice(field, inputs[1], glvBasis)
		sp := ecc.SplitScalar(inputs[0], glvBasis)
		outputs[0].Set(&(sp[0]))
		outputs[1].Set(&(sp[1]))
		// we need the absolute values for the in-circuit computations,
		// otherwise the negative values will be reduced modulo the SNARK scalar
		// field and not the emulated field.
		// 		output0 = |s0| mod r
		// 		output1 = |s1| mod r
		if outputs[0].Sign() == -1 {
			outputs[0].Neg(outputs[0])
		}
		if outputs[1].Sign() == -1 {
			outputs[1].Neg(outputs[1])
		}

		return nil
	})
}

func decomposeScalarG1Signs(mod *big.Int, inputs []*big.Int, outputs []*big.Int) error {
	return emulated.UnwrapHintWithNativeOutput(inputs, outputs, func(field *big.Int, inputs, outputs []*big.Int) error {
		if len(inputs) != 2 {
			return errors.New("expecting two inputs")
		}
		if len(outputs) != 2 {
			return errors.New("expecting two outputs")
		}
		glvBasis := new(ecc.Lattice)
		ecc.PrecomputeLattice(field, inputs[1], glvBasis)
		sp := ecc.SplitScalar(inputs[0], glvBasis)
		outputs[0].SetUint64(0)
		if sp[0].Sign() == -1 {
			outputs[0].SetUint64(1)
		}
		outputs[1].SetUint64(0)
		if sp[1].Sign() == -1 {
			outputs[1].SetUint64(1)
		}

		return nil
	})
}

// it currently supports only:
// BN254, BLS12-381, BW6-761 and Secp256k1, P256, P384 and STARK curve.
func scalarMulHint(_ *big.Int, inputs []*big.Int, outputs []*big.Int) error {
	return emulated.UnwrapHintWithNativeInput(inputs, outputs, func(field *big.Int, inputs, outputs []*big.Int) error {
		if len(outputs) != 2 {
			return errors.New("expecting two outputs")
		}
		if field.Cmp(elliptic.P256().Params().P) == 0 {
			var fp emparams.P256Fp
			var fr emparams.P256Fr
			PXLimbs := inputs[:fp.NbLimbs()]
			PYLimbs := inputs[fp.NbLimbs() : 2*fp.NbLimbs()]
			SLimbs := inputs[2*fp.NbLimbs():]
			Px, Py, S := new(big.Int), new(big.Int), new(big.Int)
			if err := Recompose(PXLimbs, fp.BitsPerLimb(), Px); err != nil {
				return err

			}
			if err := Recompose(PYLimbs, fp.BitsPerLimb(), Py); err != nil {
				return err

			}
			if err := Recompose(SLimbs, fr.BitsPerLimb(), S); err != nil {
				return err

			}
			curve := elliptic.P256()
			// compute the resulting point [s]P
			Qx, Qy := curve.ScalarMult(Px, Py, S.Bytes())
			outputs[0].Set(Qx)
			outputs[1].Set(Qy)
			return nil
		}
		if field.Cmp(bn_fp.Modulus()) == 0 {
			var fp emparams.BN254Fp
			var fr emparams.BN254Fr
			PXLimbs := inputs[:fp.NbLimbs()]
			PYLimbs := inputs[fp.NbLimbs() : 2*fp.NbLimbs()]
			SLimbs := inputs[2*fp.NbLimbs():]
			Px, Py, S := new(big.Int), new(big.Int), new(big.Int)
			if err := Recompose(PXLimbs, fp.BitsPerLimb(), Px); err != nil {
				return err

			}
			if err := Recompose(PYLimbs, fp.BitsPerLimb(), Py); err != nil {
				return err

			}
			if err := Recompose(SLimbs, fr.BitsPerLimb(), S); err != nil {
				return err

			}
			// compute the resulting point [s]Q
			var P bn254.G1Affine
			P.X.SetBigInt(Px)
			P.Y.SetBigInt(Py)
			P.ScalarMultiplication(&P, S)
			P.X.BigInt(outputs[0])
			P.Y.BigInt(outputs[1])
			return nil
		}
		return errors.New("unsupported curve")
	})
}

func halfGCDSigns(mod *big.Int, inputs []*big.Int, outputs []*big.Int) error {
	return emulated.UnwrapHintWithNativeOutput(inputs, outputs, func(field *big.Int, inputs, outputs []*big.Int) error {
		if len(inputs) != 1 {
			return fmt.Errorf("expecting one input")
		}
		if len(outputs) != 1 {
			return fmt.Errorf("expecting one output")
		}
		glvBasis := new(ecc.Lattice)
		ecc.PrecomputeLattice(field, inputs[0], glvBasis)
		outputs[0].SetUint64(0)
		if glvBasis.V1[1].Sign() == -1 {
			outputs[0].SetUint64(1)
		}

		return nil
	})
}

func halfGCD(mod *big.Int, inputs, outputs []*big.Int) error {
	return emulated.UnwrapHint(inputs, outputs, func(field *big.Int, inputs, outputs []*big.Int) error {
		if len(inputs) != 1 {
			return fmt.Errorf("expecting one input")
		}
		if len(outputs) != 2 {
			return fmt.Errorf("expecting two outputs")
		}
		glvBasis := new(ecc.Lattice)
		ecc.PrecomputeLattice(field, inputs[0], glvBasis)
		outputs[0].Set(&glvBasis.V1[0])
		outputs[1].Set(&glvBasis.V1[1])

		// we need the absolute values for the in-circuit computations,
		// otherwise the negative values will be reduced modulo the SNARK scalar
		// field and not the emulated field.
		// 		output0 = |s0| mod r
		// 		output1 = |s1| mod r
		if outputs[1].Sign() == -1 {
			outputs[1].Neg(outputs[1])
		}

		return nil
	})
}
