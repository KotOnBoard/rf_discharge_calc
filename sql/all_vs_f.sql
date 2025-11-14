SELECT
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
    MAX(CASE WHEN gas = 'Ar' THEN (Г1_B_eff-Г1_Mo_eff)/(Г1_B_eff+Г1_Mo_eff) END) as s_ar,
    MAX(CASE WHEN gas = 'Ne' THEN (Г1_B_eff-Г1_Mo_eff)/(Г1_B_eff+Г1_Mo_eff) END) as s_ne,
    MAX(CASE WHEN gas = 'He' THEN (Г1_B_eff-Г1_Mo_eff)/(Г1_B_eff+Г1_Mo_eff) END) as s_he
FROM test_table 
WHERE f IS NOT NULL 
    AND pwr = 200.0 
    AND assy = 8 
    AND S1 = 0.04625 
    AND p = 2
    AND gas IN ('Ar', 'Ne', 'He')
GROUP BY f
ORDER BY f;
