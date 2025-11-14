#!/bin/bash

# Execute SQL from file and export to CSV
query2csv() {
    local db_name="$1"
    local sql_file="$2"
    local output_csv="./csv/$3"
    
    # Validate files
    if [[ ! -f "$db_name" ]]; then
        echo "Error: Database file '$db_name' not found!"
        exit 1
    fi
    
    if [[ ! -f "$sql_file" ]]; then
        echo "Error: SQL file '$sql_file' not found!"
        exit 1
    fi
    
    # Read query from file
    query=$(cat "$sql_file")
    
    # Execute and export
    sqlite3 -header -csv "$db_name" "$query" > "$output_csv"
    
    if [[ $? -eq 0 ]]; then
        echo "Query from $sql_file executed successfully!"
        echo "Results exported to $output_csv"
    else
        echo "Error: Failed to execute query!"
        exit 1
    fi
}

query2csv "../output/test.db" "all_vs_ppl.sql" "all_vs_ppl.csv"
query2csv "../output/test.db" "all_vs_f.sql" "all_vs_f.csv"
query2csv "../output/test.db" "all_vs_ppl_vs_f.sql" "all_vs_ppl_vs_f.csv"
query2csv "../output/test.db" "vb_vmo.sql" "nc.csv"

gnuplot "ub_urf_vs_ppl.gp"
gnuplot "ub_urf_vs_f.gp"
gnuplot "d1h_b_s_vs_f.gp"
gnuplot "d1h_b_s_vs_ppl.gp"
gnuplot "ub_vs_ppl_vs_f.gp"
gnuplot "all_vs_f.gp"
gnuplot "all_vs_ppl.gp"
gnuplot "nc.gp"
gnuplot "nc2.gp"
gnuplot "nc3.gp"
