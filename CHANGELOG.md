# CHANGELOG


## v5.0.0 (2025-03-26)

### Bug Fixes

- Remove variant_record from to_dataframe()
  ([`015f823`](https://github.com/EUCANCan/variant-extractor/commit/015f823f4960a9e6af9fa6c34dcbcf4dec125d93))

BREAKING CHANGE: variant_record is no longer returned from to_dataframe() to reduce memory footprint

### Breaking Changes

- Variant_record is no longer returned from to_dataframe() to reduce memory footprint


## v4.0.10 (2025-03-26)

### Bug Fixes

- Add pysam missing dependency
  ([`0660910`](https://github.com/EUCANCan/variant-extractor/commit/0660910e5ca16e2df44d5c8cba4b905430aeabb3))

- Reduce dataframe memory footprint
  ([`d81bcd9`](https://github.com/EUCANCan/variant-extractor/commit/d81bcd90ee390bbb618f858350fc949ade131a27))


## v4.0.9 (2025-03-26)

### Bug Fixes

- Allow for numbers in shorthand SVs
  ([`8feb749`](https://github.com/EUCANCan/variant-extractor/commit/8feb7494587b0517d26c63f5e25d480d148a3771))


## v4.0.8 (2024-07-17)

### Bug Fixes

- Dummy release to generate Zenodo DOI
  ([`2c890f9`](https://github.com/EUCANCan/variant-extractor/commit/2c890f99f006f74d247262d6cd3e7995e8bb2b16))


## v4.0.7 (2024-06-28)

### Bug Fixes

- Remove author from __init__.py
  ([`2d48c4a`](https://github.com/EUCANCan/variant-extractor/commit/2d48c4a3bcdabb6e2d1842b562ac25f0d2025faa))

### Performance Improvements

- Reduce dataframe memory footprint
  ([`011f2bb`](https://github.com/EUCANCan/variant-extractor/commit/011f2bb2612f1b2f7f38701e95919f2b63c3f222))


## v4.0.6 (2023-03-08)

### Bug Fixes

- Fix missing fields in str conversion
  ([`4bbef95`](https://github.com/EUCANCan/variant-extractor/commit/4bbef9546c1c3992a75f3146f67dac0def76c641))


## v4.0.5 (2023-03-06)

### Bug Fixes

- Fix Python3.6 typings
  ([`5b3bdd6`](https://github.com/EUCANCan/variant-extractor/commit/5b3bdd696437747c5a9cf3874392817bdee0f2d7))


## v4.0.4 (2023-03-02)

### Performance Improvements

- Improve __str__ performance
  ([`ec14074`](https://github.com/EUCANCan/variant-extractor/commit/ec140748b5abd31284580c5eacfd7bb2e44cacde))


## v4.0.3 (2023-03-02)

### Performance Improvements

- Improve I/O speed by lazy loading fields
  ([`93782fc`](https://github.com/EUCANCan/variant-extractor/commit/93782fc05585981a3d3ef20a6705e18160486a47))


## v4.0.2 (2023-02-16)

### Bug Fixes

- Add __init__ file to private module
  ([`549b7c5`](https://github.com/EUCANCan/variant-extractor/commit/549b7c53fcd6c8ede96c177c4f65e414b260c7f5))


## v4.0.1 (2023-02-15)

### Bug Fixes

- Add PyPI docs
  ([`fd3b89c`](https://github.com/EUCANCan/variant-extractor/commit/fd3b89c995810834c13389f1ee16d0b9333f43ba))


## v4.0.0 (2023-02-10)

### Bug Fixes

- Change SNV length to 1
  ([`dfe6463`](https://github.com/EUCANCan/variant-extractor/commit/dfe6463e018a1d0368668e091be26a27f4dcdda7))

- Fix float .2f VCF representation
  ([`ed34a40`](https://github.com/EUCANCan/variant-extractor/commit/ed34a40bafce79c49be4b46fe9c032478106233b))

- Fix typing for Python3.6
  ([`1648bb7`](https://github.com/EUCANCan/variant-extractor/commit/1648bb744ed911b72d421ad9340570493d8baec8))

- Fix typing Python3.6
  ([`1d72f80`](https://github.com/EUCANCan/variant-extractor/commit/1d72f80f54ffb3e6db4769f4e2f3388ba96f825b))

- Match TRA to VCF specification
  ([`1f7dd1e`](https://github.com/EUCANCan/variant-extractor/commit/1f7dd1eb99f425e787a6d0cfb75c1641316c5b8c))

BREAKING CHANGE: Translocations are no longer named TRN, but TRA instead

- Print format values in the correct order
  ([`ee3d399`](https://github.com/EUCANCan/variant-extractor/commit/ee3d399a68601d37da9391e91261feb2d0cffc8d))

- Prioritize MATEID over coords for matching breakends
  ([`c6d21bf`](https://github.com/EUCANCan/variant-extractor/commit/c6d21bf51246315adca68ed4763c831742d844a2))

- Rename to breakend notation
  ([`2d86782`](https://github.com/EUCANCan/variant-extractor/commit/2d86782a4870376d6ab2985ff477ce72eed6dd96))

- Trim extra \t
  ([`ef55c36`](https://github.com/EUCANCan/variant-extractor/commit/ef55c36f05e3c8c3793a5d8b6989c5d304a640f2))

- Unify different filters for paired breakends
  ([`2407216`](https://github.com/EUCANCan/variant-extractor/commit/24072160ac4f05d8d89928df00c4fa00125bf766))

### Features

- Add empty_dataframe method
  ([`ebaf357`](https://github.com/EUCANCan/variant-extractor/commit/ebaf357f40f3dd1ac2e5bbd6c74615385e668e12))

- Add to_dataframe() method
  ([`4e666aa`](https://github.com/EUCANCan/variant-extractor/commit/4e666aa6821f58b887feb8f013afc76060d1aaff))

- Convert 1-length DELs to INS
  ([`ae6cd8c`](https://github.com/EUCANCan/variant-extractor/commit/ae6cd8ca66205a40e5e39bef06ede0949a20a0aa))

### Breaking Changes

- Translocations are no longer named TRN, but TRA instead


## v3.1.0 (2022-06-14)

### Bug Fixes

- Add missing END to info str
  ([`3606000`](https://github.com/EUCANCan/variant-extractor/commit/36060003088c8240c4d02c733a124785822a922b))

- Correctly display flag INFO values
  ([`bd83771`](https://github.com/EUCANCan/variant-extractor/commit/bd8377159d70eb6c75b181b0c35a04bd4508d102))

- Fix permutation of SV variant
  ([`d1cad92`](https://github.com/EUCANCan/variant-extractor/commit/d1cad92800dc6e603a90800943449180700b2d09))

- Remove 'None' from sample field
  ([`23bc5f3`](https://github.com/EUCANCan/variant-extractor/commit/23bc5f3caca4fc5114299480362f23bc6284ecb5))

- Remove complex indel normalization
  ([`65c5007`](https://github.com/EUCANCan/variant-extractor/commit/65c5007f4c5314abeb3baadd899ac0fd334f2740))

- Set dot for None values in sample field
  ([`43e4c83`](https://github.com/EUCANCan/variant-extractor/commit/43e4c831f5eaddb5ed2233b3d544070941be805d))

### Features

- Replace N with correct nucleotide if possible
  ([`8a260e1`](https://github.com/EUCANCan/variant-extractor/commit/8a260e10bc0f01f268a27dba3e4ed85adbc0f042))


## v3.0.0 (2022-04-25)

### Bug Fixes

- Add dot if info is not defined
  ([`70fa09b`](https://github.com/EUCANCan/variant-extractor/commit/70fa09b7d7b9eef63a826db0f1749fc850a04e8d))

### Features

- Add __str__ method to VariantRecord
  ([`7f3aa7d`](https://github.com/EUCANCan/variant-extractor/commit/7f3aa7d312c99c3654ac22b071c2468cb7432c6b))

- Read_vcf() returns a lazy generator
  ([`cd7decf`](https://github.com/EUCANCan/variant-extractor/commit/cd7decf68e18f0bf73568db7ac3e5b0b317ef65d))

BREAKING CHANGE: VariantExtractor constructor parameter only_pass renamed to pass_only

BREAKING CHANGE: read_vcf() no longer returns a list, but a lazy generator

- Switch to VariantExtractor iterations through __iter__
  ([`b7d8e58`](https://github.com/EUCANCan/variant-extractor/commit/b7d8e5846979c696dd587ec2e747f782d0c8360c))

BREAKING CHANGE: read_vcf() method no longer exists. Iterating through the VariantExtractor instance
  have the same effect.

### Breaking Changes

- Read_vcf() method no longer exists. Iterating through the VariantExtractor instance have the same
  effect.


## v2.1.0 (2022-04-20)

### Bug Fixes

- Fix imprecise breakends handling
  ([`b2cab22`](https://github.com/EUCANCan/variant-extractor/commit/b2cab22f1dc476c9e703537f389a66b94f5a7a4c))

- Fix mate_id referenced before assignment
  ([`5efdeca`](https://github.com/EUCANCan/variant-extractor/commit/5efdeca18f742a4c6453eb91916954ec600a78e6))

- Fix multiallelic record handling
  ([`0e1620f`](https://github.com/EUCANCan/variant-extractor/commit/0e1620f4dd4d08bda18399418aa7e6dcc26feaee))

- Svlen not set now defaults to 0
  ([`6afbb0d`](https://github.com/EUCANCan/variant-extractor/commit/6afbb0df1410859043e5645b80680763cf0f5310))

- Switch regex to fullmatch
  ([`c554b52`](https://github.com/EUCANCan/variant-extractor/commit/c554b5247621d95a4a81470683ee668311bdfea3))

### Features

- Add support for multiallelic variants
  ([`c763b97`](https://github.com/EUCANCan/variant-extractor/commit/c763b97b613e9cb85a98401e09888f4d5322040e))

- Divide compound indels
  ([`96f72e5`](https://github.com/EUCANCan/variant-extractor/commit/96f72e51dbea136b2f1c638b6d686a57a294e187))

- Include all information in VariantRecord
  ([`894763e`](https://github.com/EUCANCan/variant-extractor/commit/894763e38b5ca34abe3d909fe1cc2e738133d4f0))

- Infer TRN brackets from pair
  ([`6cd0ede`](https://github.com/EUCANCan/variant-extractor/commit/6cd0edef719d8c817b071331c9a3af598925a3e4))


## v2.0.0 (2022-04-04)

### Bug Fixes

- Add support for dots in bracket notation
  ([`e4dc00b`](https://github.com/EUCANCan/variant-extractor/commit/e4dc00bbd634ff7e9f343a43c365202ac02dbcc1))

- Complete the documentation
  ([`fbe5c4d`](https://github.com/EUCANCan/variant-extractor/commit/fbe5c4dbfe7df1c9fd17a06637cb76e086236add))

- Divide INV into bracket notation
  ([`c18a40c`](https://github.com/EUCANCan/variant-extractor/commit/c18a40c3af15b9935150fedbe62b8bad40aa72d7))

- Fix import path
  ([`a9db147`](https://github.com/EUCANCan/variant-extractor/commit/a9db147fe778cca3c8bc8e7a7e68cff9063ccd4a))

- Fix SGL regex
  ([`b58c867`](https://github.com/EUCANCan/variant-extractor/commit/b58c867597a86bddee0df72e6ee18baa1868b114))

### Features

- Add `qual` field
  ([`4d23686`](https://github.com/EUCANCan/variant-extractor/commit/4d236869d1457c504860b9e5b482215fc45475be))

- Add length field to VariantRecord
  ([`9a0e446`](https://github.com/EUCANCan/variant-extractor/commit/9a0e446a307c880221b39849b91c0c9db42b9c25))

- Add only_pass parameter
  ([`39d7b41`](https://github.com/EUCANCan/variant-extractor/commit/39d7b4120ce9dd389febc273a91f970af5b8bf30))

- Add support for SGLs ([#1](https://github.com/EUCANCan/variant-extractor/pull/1),
  [`f3621c0`](https://github.com/EUCANCan/variant-extractor/commit/f3621c0c9a6b6ea673666d61094f157a09b1acbe))

- Change alts to alt
  ([`8e1dd86`](https://github.com/EUCANCan/variant-extractor/commit/8e1dd86e8bae8094e047278a6c8fe60efb46eba6))

BREAKING CHANGE: VariantRecord.alts no longer contains a list, it has been changed to
  VariantRecord.alt

- Import implementation
  ([`8a676e7`](https://github.com/EUCANCan/variant-extractor/commit/8a676e78939c555c157f6090cafa74581d3b1d52))

- Include variant_type in VariantRecord
  ([`2fe8b37`](https://github.com/EUCANCan/variant-extractor/commit/2fe8b3702f1ababdb6c75290651b1f3f48f0e0a4))

BREAKING CHANGE: read_vcf() no longer returns a tuple of VariantType and VariantRecord. It now
  returns a list of just VariantRecord

- Remove INDEL_DEL and INDEL_INS
  ([`e7bbd40`](https://github.com/EUCANCan/variant-extractor/commit/e7bbd40afc84178e02372957b4060ec853e1a297))

BREAKING CHANGE: Now INDEL_DEL and INDEL_INS are treated like DEL and INS
