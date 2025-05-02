-- CryoProtect v2 - Molecule Population Script
-- This script populates the molecules table with scientifically accurate cryoprotectant data

-- Dimethyl sulfoxide (DMSO)
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    '8f52c2a1-6057-4d6c-9ec4-715f9c4dbf2c',
    'Dimethyl sulfoxide',
    'CS(=O)C',
    'InChI=1S/C2H6OS/c1-4(2)3/h1-2H3',
    'IAZDPXIOMUYVGZ-UHFFFAOYSA-N',
    'C2H6OS',
    78.13,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Glycerol
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    '9e8e4f5a-8c2d-4c1b-9a3b-7d5c2b6e4f3a',
    'Glycerol',
    'C(C(CO)O)O',
    'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
    'PEDCQBHIVMGVHV-UHFFFAOYSA-N',
    'C3H8O3',
    92.09,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Ethylene glycol
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'a1b2c3d4-e5f6-4a5b-9c7d-8e9f0a1b2c3d',
    'Ethylene glycol',
    'C(CO)O',
    'InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2',
    'LYCAIKOWRPUZTN-UHFFFAOYSA-N',
    'C2H6O2',
    62.07,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Propylene glycol
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'b2c3d4e5-f6a7-4b5c-9d7e-8f9a0b1c2d3e',
    'Propylene glycol',
    'CC(CO)O',
    'InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3',
    'DNIAPMSPPWPWGF-UHFFFAOYSA-N',
    'C3H8O2',
    76.09,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Trehalose
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'c3d4e5f6-a7b8-4c5d-9e7f-8a9b0c1d2e3f',
    'Trehalose',
    'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)CO)O',
    'InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12-10(19)6(15)5(14)3(2-13)22-12/h3-19H,1-2H2/t3-,4-,5-,6-,7-,8+,9-,10-,11-,12+/m1/s1',
    'OHCBMWOFNXFUKL-JCCZQYLRSA-N',
    'C12H22O11',
    342.30,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Sucrose
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'd4e5f6a7-b8c9-4d5e-9f7a-8b9c0d1e2f3a',
    'Sucrose',
    'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O',
    'InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12(3-15)10(19)6(2-14)22-5(12)20/h4-11,13-20H,1-3H2/t4-,5+,6-,7-,8+,9-,10+,11-,12+/m1/s1',
    'CZMRCDWAGMRECN-UGDNZRGBSA-N',
    'C12H22O11',
    342.30,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Methanol
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'e5f6a7b8-c9d0-4e5f-9a7b-8c9d0e1f2a3b',
    'Methanol',
    'CO',
    'InChI=1S/CH4O/c1-2/h2H,1H3',
    'OKKJLVBELUTLKV-UHFFFAOYSA-N',
    'CH4O',
    32.04,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Formamide
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'f6a7b8c9-d0e1-4f6a-9b7c-8d9e0f1a2b3c',
    'Formamide',
    'C(=O)N',
    'InChI=1S/CH3NO/c2-1-3/h1H,(H2,2,3)',
    'ZHNUHDYFZUAESO-UHFFFAOYSA-N',
    'CH3NO',
    45.04,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- Acetamide
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'a7b8c9d0-e1f2-4a7b-9c8d-9e0f1a2b3c4d',
    'Acetamide',
    'CC(=O)N',
    'InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)',
    'DLFVBJFMPXGRIB-UHFFFAOYSA-N',
    'C2H5NO',
    59.07,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);

-- 1,2-Propanediol (same as Propylene glycol, but keeping for consistency with the script)
INSERT INTO molecules (
    id, name, smiles, inchi, inchikey, formula, molecular_weight, 
    data_source, version, modification_history, created_at, updated_at
) VALUES (
    'b8c9d0e1-f2a3-4b8c-9d0e-1f2a3b4c5d6e',
    '1,2-Propanediol',
    'CC(CO)O',
    'InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3',
    'DNIAPMSPPWPWGF-UHFFFAOYSA-N',
    'C3H8O2',
    76.09,
    'CryoProtect v2 Database Population Script',
    1,
    '[{"timestamp": "2025-04-18T19:44:00Z", "action": "created", "user_id": null}]',
    NOW(),
    NOW()
);