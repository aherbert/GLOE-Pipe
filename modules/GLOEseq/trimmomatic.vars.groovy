//vars for task Trimmomatic 
Trimmomatic_THREADS=" -threads " + Integer.toString(ESSENTIAL_THREADS) // threads to use
Trimmomatic_ADAPTER=ESSENTIAL_TRIMMOMATIC
Trimmomatic_MINREADLEN=" MINLEN:" + ESSENTIAL_MIN_READLENGTH