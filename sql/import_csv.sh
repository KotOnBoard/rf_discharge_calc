#!/bin/bash

csv_to_sqlite() {
    local csv="$1"
    local db="${2:-data.db}"
    local table="${3:-$(basename "$csv" .csv)}"
    
    echo "Importing $csv into $db as $table..."
    
    # Import and let SQLite auto-detect types
CREATE TABLE "test_table" (
	"ro_B"	TEXT,
	"M_B"	TEXT,
	"ro_Mo"	TEXT,
	"M_Mo"	TEXT,
	"ro_MoO3"	TEXT,
	"M_MoO3"	TEXT,
	"M_Si"	TEXT,
	"ro_Si"	TEXT,
	"e"	TEXT,
	"sm"	TEXT,
	"k"	TEXT,
	"e0"	TEXT,
	"N_A"	TEXT,
	"Ti"	TEXT,
	"Tn"	TEXT,
	"m"	TEXT,
	"gas"	TEXT,
	"l"	REAL,
	"f"	REAL,
	"Pwr"	REAL,
	"Assy"	REAL,
	"p"	REAL,
	"S1"	REAL,
	"R"	TEXT,
	"S2"	TEXT,
	"omega"	TEXT,
	"d"	TEXT,
	"ng"	TEXT,
	"lam_i"	TEXT,
	"lam_i_d"	TEXT,
	"hl"	TEXT,
	"hR"	TEXT,
	"deff"	TEXT,
	"Te"	TEXT,
	"ub"	TEXT,
	"Kiz"	TEXT,
	"Kel"	TEXT,
	"Kex"	TEXT,
	"vm"	TEXT,
	"xi_c"	TEXT,
	"Vf"	TEXT,
	"sm2"	TEXT,
	"sm1"	TEXT,
	"C1"	TEXT,
	"C2"	TEXT,
	"Vp_amp"	TEXT,
	"Vp_avg"	TEXT,
	"Vp_symm"	TEXT,
	"Ubias"	TEXT,
	"V1_avg"	TEXT,
	"V2_avg"	TEXT,
	"Sohm"	TEXT,
	"Sstoc1"	TEXT,
	"Sstoc2"	TEXT,
	"ns1_55"	TEXT,
	"n0"	TEXT,
	"V"	TEXT,
	"Ji_symm"	TEXT,
	"Ji1"	TEXT,
	"J1"	TEXT,
	"ns2_61"	TEXT,
	"Ji2"	TEXT,
	"Sabs"	TEXT,
	"ns"	TEXT,
	"Vrf"	TEXT,
	"dEi1"	TEXT,
	"dEi1 [eV]"	TEXT,
	"dEi2"	TEXT,
	"dEi2 [eV]"	TEXT,
	"dEi_symm"	TEXT,
	"dEi_symm [eV]"	TEXT,
	"Ns1_72"	TEXT,
	"Ns2_73"	TEXT,
	"NN1"	TEXT,
	"G1_Mo_eff"	TEXT,
	"G1_B_eff"	TEXT,
	"G1_Si_eff"	TEXT,
	"NN2"	TEXT,
	"G2_Mo_eff"	TEXT,
	"G2_B_eff"	TEXT,
	"G2_Si_eff"	TEXT,
	"d1_MoB"	TEXT,
	"d1M_Mo"	TEXT,
	"d1H_Mo"	TEXT,
	"d1M_B"	TEXT,
	"d1H_B"	TEXT,
	"s1_B"	TEXT,
	"d1M_Si"	TEXT,
	"d1H_Si"	TEXT,
	"s1_Si"	TEXT,
	"d2_MoB"	TEXT,
	"d2M_Mo"	TEXT,
	"d2H_Mo"	TEXT,
	"d2M_B"	TEXT,
	"d2H_B"	TEXT,
	"s2_B"	TEXT,
	"d2M_Si"	TEXT,
	"d2H_Si"	TEXT,
	"s2_Si"	TEXT
);
    sqlite3 "$db" ".mode csv" ".import '$csv' '$table'"
    
    # The import automatically handles:
    # - Integer detection (whole numbers)
    # - Real detection (numbers with decimals)  
    # - Text detection (everything else)
    
    echo "Done! Schema:"
    sqlite3 "$db" ".schema $table"
}

csv_to_sqlite "../output/HighDef+f.csv" "LowDef.db" "test_table"
