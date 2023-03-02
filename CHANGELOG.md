# Changelog

<!--next-version-placeholder-->

## v4.0.4 (2023-03-02)
### Performance
* Improve __str__ performance ([`ec14074`](https://github.com/EUCANCan/variant-extractor/commit/ec140748b5abd31284580c5eacfd7bb2e44cacde))

## v4.0.3 (2023-03-02)
### Performance
* Improve I/O speed by lazy loading fields ([`93782fc`](https://github.com/EUCANCan/variant-extractor/commit/93782fc05585981a3d3ef20a6705e18160486a47))

## v4.0.2 (2023-02-16)
### Fix
* Add __init__ file to private module ([`549b7c5`](https://github.com/EUCANCan/variant-extractor/commit/549b7c53fcd6c8ede96c177c4f65e414b260c7f5))

## v4.0.1 (2023-02-15)
### Fix
* Add PyPI docs ([`fd3b89c`](https://github.com/EUCANCan/variant-extractor/commit/fd3b89c995810834c13389f1ee16d0b9333f43ba))

## v4.0.0 (2023-02-10)
### Feature
* Add empty_dataframe method ([`ebaf357`](https://github.com/EUCANCan/variant-extractor/commit/ebaf357f40f3dd1ac2e5bbd6c74615385e668e12))
* Add to_dataframe() method ([`4e666aa`](https://github.com/EUCANCan/variant-extractor/commit/4e666aa6821f58b887feb8f013afc76060d1aaff))
* Convert 1-length DELs to INS ([`ae6cd8c`](https://github.com/EUCANCan/variant-extractor/commit/ae6cd8ca66205a40e5e39bef06ede0949a20a0aa))

### Fix
* Match TRA to VCF specification ([`1f7dd1e`](https://github.com/EUCANCan/variant-extractor/commit/1f7dd1eb99f425e787a6d0cfb75c1641316c5b8c))
* Fix typing Python3.6 ([`1d72f80`](https://github.com/EUCANCan/variant-extractor/commit/1d72f80f54ffb3e6db4769f4e2f3388ba96f825b))
* Fix typing for Python3.6 ([`1648bb7`](https://github.com/EUCANCan/variant-extractor/commit/1648bb744ed911b72d421ad9340570493d8baec8))
* Rename to breakend notation ([`2d86782`](https://github.com/EUCANCan/variant-extractor/commit/2d86782a4870376d6ab2985ff477ce72eed6dd96))
* Print format values in the correct order ([`ee3d399`](https://github.com/EUCANCan/variant-extractor/commit/ee3d399a68601d37da9391e91261feb2d0cffc8d))
* Fix wrong SAMPLEs in multiallelic variants ([`7fc170f`](https://github.com/EUCANCan/variant-extractor/commit/7fc170fd535b2483f8d344428a5e567148e8b4c6))
* Fix float .2f VCF representation ([`ed34a40`](https://github.com/EUCANCan/variant-extractor/commit/ed34a40bafce79c49be4b46fe9c032478106233b))
* Trim extra \t ([`ef55c36`](https://github.com/EUCANCan/variant-extractor/commit/ef55c36f05e3c8c3793a5d8b6989c5d304a640f2))
* Prioritize MATEID over coords for matching breakends ([`c6d21bf`](https://github.com/EUCANCan/variant-extractor/commit/c6d21bf51246315adca68ed4763c831742d844a2))
* Unify different filters for paired breakends ([`2407216`](https://github.com/EUCANCan/variant-extractor/commit/24072160ac4f05d8d89928df00c4fa00125bf766))
* Change SNV length to 1 ([`dfe6463`](https://github.com/EUCANCan/variant-extractor/commit/dfe6463e018a1d0368668e091be26a27f4dcdda7))

### Breaking
* Translocations are no longer named TRN, but TRA instead ([`1f7dd1e`](https://github.com/EUCANCan/variant-extractor/commit/1f7dd1eb99f425e787a6d0cfb75c1641316c5b8c))

## v3.1.0 (2022-06-14)
### Feature
* Replace N with correct nucleotide if possible ([`8a260e1`](https://github.com/EUCANCan/variant-extractor/commit/8a260e10bc0f01f268a27dba3e4ed85adbc0f042))

### Fix
* Remove 'None' from sample field ([`23bc5f3`](https://github.com/EUCANCan/variant-extractor/commit/23bc5f3caca4fc5114299480362f23bc6284ecb5))
* Correctly display flag INFO values ([`bd83771`](https://github.com/EUCANCan/variant-extractor/commit/bd8377159d70eb6c75b181b0c35a04bd4508d102))
* Fix permutation of SV variant ([`d1cad92`](https://github.com/EUCANCan/variant-extractor/commit/d1cad92800dc6e603a90800943449180700b2d09))
* Add missing END to info str ([`3606000`](https://github.com/EUCANCan/variant-extractor/commit/36060003088c8240c4d02c733a124785822a922b))
* Remove complex indel normalization ([`65c5007`](https://github.com/EUCANCan/variant-extractor/commit/65c5007f4c5314abeb3baadd899ac0fd334f2740))
* Set dot for None values in sample field ([`43e4c83`](https://github.com/EUCANCan/variant-extractor/commit/43e4c831f5eaddb5ed2233b3d544070941be805d))

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
