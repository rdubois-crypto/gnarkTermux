# Example

We create a FakeGLV circuit, where $[k]P = Q$ is verified by decomposing $k$ into a fraction of small sub-scalars.

The test is done  with P256 (GLV is not available for this curve).


```
go test -v -run TestScalarMulFakeGLV
```
Output:
```
=== RUN   TestScalarMulFakeGLV
Circuit is compiled.
Setup is done (6.0).
Proof is generated (0.3).
Verification succeeded (0.0)!
--- PASS: TestScalarMulFakeGLV (7.19s)
PASS
ok  	github.com/simonmasson/GLVInCircuit/go/circuits	7.236s
```