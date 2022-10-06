# curdleproofs.pie: Verifier Only version

curdleproofs.pie is a Python implementation of the [curdleproofs protocol](https://github.com/asn-d6/curdleproofs) from the Ethereum Foundation Cryptography Research team.

The implementation is at feature parity with the original implementation and passes all the tests provided there.

Additionally, curdleproofs.pie contains a self-contained package implementing [Merlin Transcripts](https://merlin.cool) (along with the required subset of the STROBE framework) in Python. This is also tested for equivalence to the [Rust implementation](https://crates.io/crates/merlin) using the provided tests in that repository.

This is a branch where the prover code is stripped and all that is left is serde and verification code. This code is relevant to consensus clients for verifying and replicating in clients.

See an example of the verifier usage in [`test_curdleproofs.py`](https://github.com/nalinbhardwaj/curdleproofpie/blob/verifier-only/curdleproofs/test_curdleproofs.py).

Note that the [`master`](https://github.com/nalinbhardwaj/curdleproofs.pie) branch additionally provides a prover implementation compatible with this verifier.
