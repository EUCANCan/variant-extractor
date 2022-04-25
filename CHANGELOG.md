# Changelog

<!--next-version-placeholder-->

## v3.0.0 (2022-04-25)
### Feature
* Switch to VariantExtractor iterations through __iter__ ([`b7d8e58`](https://github.com/EUCANCan/variant-extractor/commit/b7d8e5846979c696dd587ec2e747f782d0c8360c))
* Read_vcf() returns a lazy generator ([`cd7decf`](https://github.com/EUCANCan/variant-extractor/commit/cd7decf68e18f0bf73568db7ac3e5b0b317ef65d))
* Add __str__ method to VariantRecord ([`7f3aa7d`](https://github.com/EUCANCan/variant-extractor/commit/7f3aa7d312c99c3654ac22b071c2468cb7432c6b))

### Fix
* Add dot if info is not defined ([`70fa09b`](https://github.com/EUCANCan/variant-extractor/commit/70fa09b7d7b9eef63a826db0f1749fc850a04e8d))

### Breaking
* read_vcf() method no longer exists. Iterating through the VariantExtractor instance have the same effect.  ([`b7d8e58`](https://github.com/EUCANCan/variant-extractor/commit/b7d8e5846979c696dd587ec2e747f782d0c8360c))
* VariantExtractor constructor parameter only_pass renamed to pass_only ([`cd7decf`](https://github.com/EUCANCan/variant-extractor/commit/cd7decf68e18f0bf73568db7ac3e5b0b317ef65d))
* read_vcf() no longer returns a list, but a lazy generator  ([`cd7decf`](https://github.com/EUCANCan/variant-extractor/commit/cd7decf68e18f0bf73568db7ac3e5b0b317ef65d))

## v2.1.0 (2022-04-20)
### Feature
* Include all information in VariantRecord ([`894763e`](https://github.com/EUCANCan/variant-extractor/commit/894763e38b5ca34abe3d909fe1cc2e738133d4f0))
* Infer TRN brackets from pair ([`6cd0ede`](https://github.com/EUCANCan/variant-extractor/commit/6cd0edef719d8c817b071331c9a3af598925a3e4))
* Add support for multiallelic variants ([`c763b97`](https://github.com/EUCANCan/variant-extractor/commit/c763b97b613e9cb85a98401e09888f4d5322040e))
* Divide compound indels ([`96f72e5`](https://github.com/EUCANCan/variant-extractor/commit/96f72e51dbea136b2f1c638b6d686a57a294e187))

### Fix
* Fix imprecise breakends handling ([`b2cab22`](https://github.com/EUCANCan/variant-extractor/commit/b2cab22f1dc476c9e703537f389a66b94f5a7a4c))
* Switch regex to fullmatch ([`c554b52`](https://github.com/EUCANCan/variant-extractor/commit/c554b5247621d95a4a81470683ee668311bdfea3))
* Fix multiallelic record handling ([`0e1620f`](https://github.com/EUCANCan/variant-extractor/commit/0e1620f4dd4d08bda18399418aa7e6dcc26feaee))
* SVLEN not set now defaults to 0 ([`6afbb0d`](https://github.com/EUCANCan/variant-extractor/commit/6afbb0df1410859043e5645b80680763cf0f5310))
* Fix mate_id referenced before assignment ([`5efdeca`](https://github.com/EUCANCan/variant-extractor/commit/5efdeca18f742a4c6453eb91916954ec600a78e6))

### Documentation
* Fix paragraph headers ([`7dc9e9d`](https://github.com/EUCANCan/variant-extractor/commit/7dc9e9d30bd8df648ea79be2f1aaf96ecf7e1d7e))

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
