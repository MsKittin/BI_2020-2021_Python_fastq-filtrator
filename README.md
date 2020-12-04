# BI_2020-2021_Python_fastq-filtrator

**Program for filtering .fastq files**

**Authors:** Ekaterina Kirillova (parsing)
         Maksim Serdakov (filtering)
         
**Supported arguments:**
Optional:
         --min_length INT
         Minimum length of reads to pass the filtering.
         --keep_filtered
         If specified - failed-filtering reads will be written to a separate file.
         --gc_bounds INT (INT)
         Takes from 1 to 2 values:
         if 1 value is specified - this is the lower bound of GC-content in read to pass the filtering,
         if 2 values are specified - these are the lower and upper bounds respectively in read to pass the filtering.
         --output_base_name STR
         Common prefix for output file(s),
         if not specified - output file(s) will be named as source file without extension.
Required:
         .fastq file
         Positional argument, must be specified last.
