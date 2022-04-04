# Changelog

<!--next-version-placeholder-->

## v2.0.0 (2022-04-04)
### Feature
* Add length field to VariantRecord ([`9a0e446`](https://github.com/EUCANCan/variant-extractor/commit/9a0e446a307c880221b39849b91c0c9db42b9c25))
* Change alts to alt ([`8e1dd86`](https://github.com/EUCANCan/variant-extractor/commit/8e1dd86e8bae8094e047278a6c8fe60efb46eba6))
* Include variant_type in VariantRecord ([`2fe8b37`](https://github.com/EUCANCan/variant-extractor/commit/2fe8b3702f1ababdb6c75290651b1f3f48f0e0a4))
* Add `qual` field ([`4d23686`](https://github.com/EUCANCan/variant-extractor/commit/4d236869d1457c504860b9e5b482215fc45475be))
* Remove INDEL_DEL and INDEL_INS ([`e7bbd40`](https://github.com/EUCANCan/variant-extractor/commit/e7bbd40afc84178e02372957b4060ec853e1a297))
* Add only_pass parameter ([`39d7b41`](https://github.com/EUCANCan/variant-extractor/commit/39d7b4120ce9dd389febc273a91f970af5b8bf30))
* Add support for SGLs ([#1](https://github.com/EUCANCan/variant-extractor/issues/1)) ([`f3621c0`](https://github.com/EUCANCan/variant-extractor/commit/f3621c0c9a6b6ea673666d61094f157a09b1acbe))

### Fix
* Divide INV into bracket notation ([`c18a40c`](https://github.com/EUCANCan/variant-extractor/commit/c18a40c3af15b9935150fedbe62b8bad40aa72d7))
* Fix SGL regex ([`b58c867`](https://github.com/EUCANCan/variant-extractor/commit/b58c867597a86bddee0df72e6ee18baa1868b114))

### Breaking
* VariantRecord.alts no longer contains a list, it has been changed to VariantRecord.alt  ([`8e1dd86`](https://github.com/EUCANCan/variant-extractor/commit/8e1dd86e8bae8094e047278a6c8fe60efb46eba6))
* read_vcf() no longer returns a tuple of VariantType and VariantRecord. It now returns a list of just VariantRecord  ([`2fe8b37`](https://github.com/EUCANCan/variant-extractor/commit/2fe8b3702f1ababdb6c75290651b1f3f48f0e0a4))
* Now INDEL_DEL and INDEL_INS are treated like DEL and INS  ([`e7bbd40`](https://github.com/EUCANCan/variant-extractor/commit/e7bbd40afc84178e02372957b4060ec853e1a297))

### Documentation
* Update documentation link ([`d488f34`](https://github.com/EUCANCan/variant-extractor/commit/d488f346f4ac171d7b8a54c85a3d54f58d6014ba))
* Update docs about imprecise breakends ([`b5b7e4d`](https://github.com/EUCANCan/variant-extractor/commit/b5b7e4db8f5f7bb030cfeade323b1e7c2b90f00e))
* Update docs with pip installation ([`d3a4a53`](https://github.com/EUCANCan/variant-extractor/commit/d3a4a5351ad9735a3feb6382ade4073a1f91a6a0))

## v1.1.0 (2022-03-23)
### Feature
* Import implementation ([`9241bc1`](https://github.com/Rapsssito/variant-extractor/commit/9241bc18298f0783718b47d10d528671ecedf30b))

### Fix
* Complete the documentation ([`a2ec02b`](https://github.com/Rapsssito/variant-extractor/commit/a2ec02b9b5af3dbbf14661841a8377bf99617b07))
* Add support for dots in bracket notation ([`37f1232`](https://github.com/Rapsssito/variant-extractor/commit/37f1232a0832e35410d0c5250a9dfdc41a361b88))
* Fix import path ([`f1c5b99`](https://github.com/Rapsssito/variant-extractor/commit/f1c5b99e92e2af2c7263cd9d89c2ac42f97aeefd))

### Documentation
* Add documentation ([`0ae6e84`](https://github.com/Rapsssito/variant-extractor/commit/0ae6e84d8d6b52bd27b44d4b1300b6148c1ce242))
* Extend docs ([`bcea586`](https://github.com/Rapsssito/variant-extractor/commit/bcea5867d2f3dc37ea7b6f6d72286e6d61d408b1))
* Add initial documentation to the README ([`5ad69bf`](https://github.com/Rapsssito/variant-extractor/commit/5ad69bfb09918fa482a46b572c0c13dd0f4a1420))

## v0.1.0 (2022-03-17)
### Feature
* Import implementation ([`9241bc1`](https://github.com/Rapsssito/variant-extractor/commit/9241bc18298f0783718b47d10d528671ecedf30b))

### Fix
* Fix import path ([`f1c5b99`](https://github.com/Rapsssito/variant-extractor/commit/f1c5b99e92e2af2c7263cd9d89c2ac42f97aeefd))

### Documentation
* Add initial documentation to the README ([`5ad69bf`](https://github.com/Rapsssito/variant-extractor/commit/5ad69bfb09918fa482a46b572c0c13dd0f4a1420))
