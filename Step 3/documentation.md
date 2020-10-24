<h1>Format</h1>
A FASTQ file normally uses four lines per sequence.

* Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description (like a FASTA title line).
* Line 2 is the raw sequence letters.
* Line 3 begins with a '+' character and is optionally followed by the same sequence identifier (and any description) again.
* Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number of symbols as letters in the sequence.
A FASTQ file containing a single sequence might look like this:

```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```
* [Wikipedia page for further reference.](https://en.wikipedia.org/wiki/FASTQ_format)
* [Read for more detailed explanation on quality value encoding](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
