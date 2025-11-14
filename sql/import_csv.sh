#!/bin/bash

csv_to_sqlite() {
    local csv="$1"
    local db="${2:-data.db}"
    local table="${3:-$(basename "$csv" .csv)}"
    
    echo "Importing $csv into $db as $table..."
    
    # Import and let SQLite auto-detect types
    sqlite3 "$db" ".mode csv" ".import '$csv' '$table'"
    
    # The import automatically handles:
    # - Integer detection (whole numbers)
    # - Real detection (numbers with decimals)  
    # - Text detection (everything else)
    
    echo "Done! Schema:"
    sqlite3 "$db" ".schema $table"
}

csv_to_sqlite "../output/TEST323CSV.csv" "test.db" "test_table"
