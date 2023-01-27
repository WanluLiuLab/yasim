import multiprocessing
import os
import shutil
import sys

from pyspark.sql import SparkSession, types

if __name__ == '__main__':
    spark = SparkSession.builder.getOrCreate()
    spark.sparkContext.setLogLevel("ERROR")
    print(f"Spark version: {spark.version}")
    FILE_DIR = os.path.dirname(os.path.abspath(__file__))

    fn = sys.argv[1]
    print(f"Reading {fn}...")
    df = spark.read.csv(
        path=fn,
        schema=types.StructType().
        add("REFERENCE_NAME", types.StringType()).
        add("REFERENCE_POS", types.IntegerType()).
        add("NUM_READS", types.IntegerType()),
        encoding="UTF-8",
        enforceSchema=True,
        inferSchema=False,
        locale="en-US",
        mode="DROPMALFORMED",
        sep="\t"
    )
    length = df.count()
    print(f"Sampling {fn} ({length} items)...")
    df = df.sample(
        withReplacement=False,
        fraction=1E-3
    ).repartition(
        multiprocessing.cpu_count()
    )
    length = df.count()
    print(f"Writing {fn} ({length} items)...")
    if os.path.exists(fn + ".sampled.parquet.d"):
        shutil.rmtree(fn + ".sampled.parquet.d")
    df.write.parquet(fn + ".sampled.parquet.d")
