package circuits

import (
	"crypto/elliptic"
	"crypto/rand"
	"fmt"
	"testing"
	"time"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/consensys/gnark/backend/groth16"
	"github.com/consensys/gnark/frontend"
	"github.com/consensys/gnark/frontend/cs/r1cs"
	"github.com/consensys/gnark/std/math/emulated"
	"github.com/consensys/gnark/test"
)

var testCurve = ecc.BN254

type ScalarMulFakeGLVTest[T, S emulated.FieldParams] struct {
	Q, R AffinePoint[T]  `gnark:",public"` // Public input
	S    emulated.Element[S]
}

func (c *ScalarMulFakeGLVTest[T, S]) Define(api frontend.API) error {
	cr, err := New[T, S](api, GetCurveParams[T]())
	if err != nil {
		return err
	}
	res := cr.scalarMulFakeGLV(&c.Q, &c.S)
	cr.AssertIsEqual(res, &c.R)
	return nil
}

func TestScalarMulFakeGLV(t *testing.T) {
	assert := test.NewAssert(t)
	p256 := elliptic.P256()
	s, err := rand.Int(rand.Reader, p256.Params().N)
	assert.NoError(err)
	px, py := p256.ScalarBaseMult(s.Bytes())

	circuit := ScalarMulFakeGLVTest[emulated.P256Fp, emulated.P256Fr]{}
	witness := ScalarMulFakeGLVTest[emulated.P256Fp, emulated.P256Fr]{
		S: emulated.ValueOf[emulated.P256Fr](s),
		Q: AffinePoint[emulated.P256Fp]{
			X: emulated.ValueOf[emulated.P256Fp](p256.Params().Gx),
			Y: emulated.ValueOf[emulated.P256Fp](p256.Params().Gy),
		},
		R: AffinePoint[emulated.P256Fp]{
			X: emulated.ValueOf[emulated.P256Fp](px),
			Y: emulated.ValueOf[emulated.P256Fp](py),
		},
	}
	err = test.IsSolved(&circuit, &witness, testCurve.ScalarField())
	assert.NoError(err)

	// Compile the circuit
	r1cs, err := frontend.Compile(ecc.BN254.ScalarField(), r1cs.NewBuilder, &circuit)
	assert.NoError(err)
	print("Circuit is compiled.\n")

	// 1. One time setup
	start := time.Now()
	pk, vk, err := groth16.Setup(r1cs)
	elapsed := time.Since(start)
	assert.NoError(err)
	fmt.Printf("Setup is done (%.1f).\n", elapsed.Seconds())

	witnessFull, err := frontend.NewWitness(&witness, ecc.BN254.ScalarField())
	assert.NoError(err)
	publicWitness, err := witnessFull.Public()
	assert.NoError(err)

	// 2. Proof creation
	start = time.Now()
	proof, err := groth16.Prove(r1cs, pk, witnessFull)
	elapsed = time.Since(start)
	assert.NoError(err)
	fmt.Printf("Proof is generated (%.1f).\n", elapsed.Seconds())

	// 3. Proof verification
	start = time.Now()
	err_verif := groth16.Verify(proof, vk, publicWitness)
	elapsed = time.Since(start)
	assert.NoError(err_verif)
	fmt.Printf("Verification succeeded (%.1f)!\n", elapsed.Seconds())

}

type ScalarMulJoyeTest[T, S emulated.FieldParams] struct {
	P, Q AffinePoint[T]
	S    emulated.Element[S]
}

func (c *ScalarMulJoyeTest[T, S]) Define(api frontend.API) error {
	cr, err := New[T, S](api, GetCurveParams[T]())
	if err != nil {
		return err
	}
	res := cr.scalarMulJoye(&c.P, &c.S)
	cr.AssertIsEqual(res, &c.Q)
	return nil
}

func TestScalarMulJoye(t *testing.T) {
	assert := test.NewAssert(t)
	p256 := elliptic.P256()
	s, err := rand.Int(rand.Reader, p256.Params().N)
	assert.NoError(err)
	px, py := p256.ScalarBaseMult(s.Bytes())

	circuit := ScalarMulJoyeTest[emulated.P256Fp, emulated.P256Fr]{}
	witness := ScalarMulJoyeTest[emulated.P256Fp, emulated.P256Fr]{
		S: emulated.ValueOf[emulated.P256Fr](s),
		P: AffinePoint[emulated.P256Fp]{
			X: emulated.ValueOf[emulated.P256Fp](p256.Params().Gx),
			Y: emulated.ValueOf[emulated.P256Fp](p256.Params().Gy),
		},
		Q: AffinePoint[emulated.P256Fp]{
			X: emulated.ValueOf[emulated.P256Fp](px),
			Y: emulated.ValueOf[emulated.P256Fp](py),
		},
	}
	err = test.IsSolved(&circuit, &witness, testCurve.ScalarField())
	assert.NoError(err)

	// Compile the circuit
	r1cs, err := frontend.Compile(ecc.BN254.ScalarField(), r1cs.NewBuilder, &circuit)
	assert.NoError(err)
	print("Circuit is compiled.\n")

	// 1. One time setup
	start := time.Now()
	pk, vk, err := groth16.Setup(r1cs)
	elapsed := time.Since(start)
	assert.NoError(err)
	fmt.Printf("Setup is done (%.1f).\n", elapsed.Seconds())

	witnessFull, err := frontend.NewWitness(&witness, ecc.BN254.ScalarField())
	assert.NoError(err)
	publicWitness, err := witnessFull.Public()
	assert.NoError(err)

	// 2. Proof creation
	start = time.Now()
	proof, err := groth16.Prove(r1cs, pk, witnessFull)
	elapsed = time.Since(start)
	assert.NoError(err)
	fmt.Printf("Proof is generated (%.1f).\n", elapsed.Seconds())

	// 3. Proof verification
	start = time.Now()
	err_verif := groth16.Verify(proof, vk, publicWitness)
	elapsed = time.Since(start)
	assert.NoError(err_verif)
	fmt.Printf("Verification succeeded (%.1f)!\n", elapsed.Seconds())

}
