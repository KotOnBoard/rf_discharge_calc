SELECT
    pwr as ppl,
	f,
    MAX(CASE WHEN gas = 'Ar' THEN Ubias END) as ub_ar,
    MAX(CASE WHEN gas = 'Ne' THEN Ubias END) as ub_ne,
    MAX(CASE WHEN gas = 'He' THEN Ubias END) as ub_he,
    MAX(CASE WHEN gas = 'Ar' THEN Vrf END) as urf_ar,
    MAX(CASE WHEN gas = 'Ne' THEN Vrf END) as urf_ne,
    MAX(CASE WHEN gas = 'He' THEN Vrf END) as urf_he,
    MAX(CASE WHEN gas = 'Ar' THEN d1H_B END) as d1h_b_ar,
    MAX(CASE WHEN gas = 'Ne' THEN d1H_B END) as d1h_b_ne,
    MAX(CASE WHEN gas = 'He' THEN d1H_B END) as d1h_b_he,
	IFNULL(MAX(CASE WHEN gas = 'Ar' THEN (d1H_B-d1H_Mo)/(d1H_B+d1H_Mo) END), -2.0) as s_ar,
	IFNULL(MAX(CASE WHEN gas = 'Ne' THEN (d1H_B-d1H_Mo)/(d1H_B+d1H_Mo) END), -2.0) as s_ne,
	IFNULL(MAX(CASE WHEN gas = 'He' THEN (d1H_B-d1H_Mo)/(d1H_B+d1H_Mo) END), -2.0) as s_he
FROM test_table 
WHERE Pwr IS NOT NULL 
    AND assy = 8 
    AND S1 = 0.04625 
    AND p = 3
    AND gas IN ('Ar', 'Ne', 'He')
GROUP BY ppl, f
ORDER BY ppl, f;