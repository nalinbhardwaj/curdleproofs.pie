# curdleproofpie: verifier only

curdleproofpie is a Python equivalent to [curdleproofs](https://github.com/asn-d6/curdleproofs) original implementation along with a python Merlin transcript implementation.

All tests pass for equivalency to the original curdleproofs Rust implementation as well as Merlin transcript conformance to Strobe.

This is a branch where the prover code is stripped and all that is left is serde and verification code. This code is relevant to consensus clients for verifying and replicating in clients.

See an example of the verifier usage in [`test_curdleproofs.py`](https://github.com/nalinbhardwaj/curdleproofpie/blob/verifier-only/curdleproofs/test_curdleproofs.py).

Note that the [master branch](https://github.com/nalinbhardwaj/curdleproofpie) additionally provides a prover implementation compatible with this verifier.
