WITH params AS (
    SELECT 50 as f_val, 2 as p_val, 150 as pwr_val, 1 as cond_name
    UNION ALL
    SELECT 70 as f_val, 2 as p_val, 150 as pwr_val, 2 as cond_name
    UNION ALL
    SELECT 50 as f_val, 2 as p_val, 250 as pwr_val, 3 as cond_name
    UNION ALL
    SELECT 70 as f_val, 2 as p_val, 250 as pwr_val, 4 as cond_name
)
SELECT 
    p.cond_name,
    MAX(CASE WHEN gas = 'Ar' THEN d1H_B/d1H_Mo END) as b2mo_ar,
    MAX(CASE WHEN gas = 'Ne' THEN d1H_B/d1H_Mo END) as b2mo_ne,
    MAX(CASE WHEN gas = 'He' THEN d1H_B/d1H_Mo END) as b2mo_he
FROM test_table t
JOIN params p ON t.f = p.f_val AND t.p = p.p_val AND t.pwr = p.pwr_val
WHERE t.assy = 8 
    AND t.gas IN ('Ar', 'Ne', 'He')
	AND S1 = 0.01
GROUP BY p.cond_name
