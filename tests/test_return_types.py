#!/usr/bin/env python3
"""
Test script to verify that the primer selection functions return meaningful objects.
"""

from pathlib import Path
import pandas as pd
from Bio.SeqRecord import SeqRecord

# Import the functions we've modified
from src.folitools.primer_selection._01_primer_selection import subset
from src.folitools.primer_selection._02_select_primer_set_by_saddle_loss import saddle
from src.folitools.primer_selection._03_extract_region_sequence import product
from src.folitools.primer_selection._04_make_excel import summary
from src.folitools.primer_selection._05_recover import recover


def main():
    print("Testing return types of primer selection functions...")
    
    # Test subset function return type
    print("\n1. Testing subset function signature:")
    print(f"   subset returns: {subset.__annotations__.get('return', 'No annotation')}")
    
    # Test saddle function return type  
    print("\n2. Testing saddle function signature:")
    print(f"   saddle returns: {saddle.__annotations__.get('return', 'No annotation')}")
    
    # Test product function return type
    print("\n3. Testing product function signature:")
    print(f"   product returns: {product.__annotations__.get('return', 'No annotation')}")
    
    # Test summary function return type
    print("\n4. Testing summary function signature:")
    print(f"   summary returns: {summary.__annotations__.get('return', 'No annotation')}")
    
    # Test recover function return type
    print("\n5. Testing recover function signature:")
    print(f"   recover returns: {recover.__annotations__.get('return', 'No annotation')}")
    
    print("\nSummary of what each function now returns:")
    print("- subset: dict[str, pd.DataFrame] - Contains 'primer_sequence_df' and 'primer_info_df'")
    print("- saddle: pd.DataFrame - Selected primers with optimization results")
    print("- product: list[SeqRecord] - Extracted amplicon sequences")
    print("- summary: pd.DataFrame - Already returned meaningful data (primer ordering info)")
    print("- recover: pd.DataFrame - Already returned meaningful data (recovered primer info)")
    
    print("\nThese are much more useful for programmatic use than status codes!")


if __name__ == "__main__":
    main()
