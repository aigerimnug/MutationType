This is an R package which takes

1. a set of mutations (single nucleotide variants, SNVs) in VCF format,
2. the corresponding reference genome (e.g., human genome hg38) and
3. a parameter “context_length” which is a positive, odd integer (e.g 3)
   
and determines for each mutation (only SNVs; other mutations like indels can be ignored) the corresponding mutation type as follows:

The mutation type is “UP[REF>ALT]DOWN” where
1. “REF>ALT” is the single nucleotide variant from REF base to ALT base, e.g., “C>T”
2. “UP” is one or more upstream bases from the reference genome (depending on the user parameter “context_length”)
3. “DOWN” is one or more downstream bases from the reference genome (same user parameter)
