# curdleproofs.pie

curdleproofs.pie is a Python implementation of the [curdleproofs protocol](https://github.com/asn-d6/curdleproofs) from the Ethereum Foundation Cryptography Research team.

The implementation is at feature parity with the original implementation and passes all the tests provided there.

Additionally, curdleproofs.pie contains a self-contained package implementing [Merlin Transcripts](https://merlin.cool) (along with the required subset of the STROBE framework) in Python. This is also tested for equivalence to the [Rust implementation](https://crates.io/crates/merlin) using the provided tests in that repository.
