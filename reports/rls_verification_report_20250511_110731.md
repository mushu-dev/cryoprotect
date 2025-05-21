# RLS Verification Report

Generated on: 2025-05-11 11:07:31

## Test Summary

- Tests run: 25
- Failures: 13
- Errors: 0
- Skipped: 3
- Success rate: 48.00%

## Failures

### test_is_project_member_invalid (tests.test_rls_policies_verification.RLSSecurityDefinerFunctionTests.test_is_project_member_invalid)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 91, in test_is_project_member_invalid
    self.assertFalse(result)
    ~~~~~~~~~~~~~~~~^^^^^^^^
AssertionError: {'error': "function is_project_member(unknown) does not exist\nLINE 1: SELECT is_project_member('22222222-2222-2222-2222-2222222222...\n               ^\nHINT:  No function matches the given name and argument types. You might need to add explicit type casts.\n"} is not false
```

### test_is_team_member (tests.test_rls_policies_verification.RLSSecurityDefinerFunctionTests.test_is_team_member)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 138, in test_is_team_member
    self.assertTrue(result)
    ~~~~~~~~~~~~~~~^^^^^^^^
AssertionError: False is not true
```

### test_user_projects (tests.test_rls_policies_verification.RLSSecurityDefinerFunctionTests.test_user_projects)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 117, in test_user_projects
    self.assertIn(self.project1_id, data_str)
    ~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: '11111111-1111-1111-1111-111111111111' not found in "[{'user_projects': '{}'}]"
```

### test_delete_molecule_as_non_member (tests.test_rls_policies_verification.RLSTableAccessPolicyTests.test_delete_molecule_as_non_member)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 451, in test_delete_molecule_as_non_member
    self.assertIn('error', result)
    ~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^
AssertionError: 'error' not found in {'count': 0, 'data': []}
```

### test_insert_molecule_as_member (tests.test_rls_policies_verification.RLSTableAccessPolicyTests.test_insert_molecule_as_member)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 291, in test_insert_molecule_as_member
    self.assertNotIn('error', result)
    ~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^
AssertionError: 'error' unexpectedly found in {'error': 'new row violates row-level security policy for table "molecules"\n'}
```

### test_select_molecule_as_member (tests.test_rls_policies_verification.RLSTableAccessPolicyTests.test_select_molecule_as_member)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 268, in test_select_molecule_as_member
    self.assertGreater(result.get('count', 0), 0)
    ~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: 0 not greater than 0
```

### test_select_molecule_as_non_member (tests.test_rls_policies_verification.RLSTableAccessPolicyTests.test_select_molecule_as_non_member)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 277, in test_select_molecule_as_non_member
    self.assertEqual(result.get('count', -1), 0)
    ~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: -1 != 0
```

### test_update_molecule_as_non_member (tests.test_rls_policies_verification.RLSTableAccessPolicyTests.test_update_molecule_as_non_member)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 380, in test_update_molecule_as_non_member
    self.assertIn('error', result)
    ~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^
AssertionError: 'error' not found in {'count': 0, 'data': []}
```

### test_experiment_property_non_user_project_access (tests.test_rls_policies_verification.RLSRelationshipPolicyTests.test_experiment_property_non_user_project_access)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 775, in test_experiment_property_non_user_project_access
    self.assertEqual(result.get('count', -1), 0)
    ~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: -1 != 0
```

### test_experiment_property_user_project_access (tests.test_rls_policies_verification.RLSRelationshipPolicyTests.test_experiment_property_user_project_access)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 766, in test_experiment_property_user_project_access
    self.assertGreater(result.get('count', 0), 0)
    ~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: 0 not greater than 0
```

### test_mixture_component_user_project_access (tests.test_rls_policies_verification.RLSRelationshipPolicyTests.test_mixture_component_user_project_access)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 749, in test_mixture_component_user_project_access
    self.assertGreater(result.get('count', 0), 0)
    ~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: 0 not greater than 0
```

### test_molecular_property_non_user_project_access (tests.test_rls_policies_verification.RLSRelationshipPolicyTests.test_molecular_property_non_user_project_access)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 741, in test_molecular_property_non_user_project_access
    self.assertEqual(result.get('count', -1), 0)
    ~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: -1 != 0
```

### test_molecular_property_user_project_access (tests.test_rls_policies_verification.RLSRelationshipPolicyTests.test_molecular_property_user_project_access)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 732, in test_molecular_property_user_project_access
    self.assertGreater(result.get('count', 0), 0)
    ~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^
AssertionError: 0 not greater than 0
```

## Skipped Tests

- test_is_project_owner (tests.test_rls_policies_verification.RLSSecurityDefinerFunctionTests.test_is_project_owner): is_project_owner function does not exist
- test_mixture_in_user_project (tests.test_rls_policies_verification.RLSSecurityDefinerFunctionTests.test_mixture_in_user_project): mixture_in_user_project function does not exist
- test_molecule_in_user_project (tests.test_rls_policies_verification.RLSSecurityDefinerFunctionTests.test_molecule_in_user_project): molecule_in_user_project function does not exist

## Log Summary

```
2025-05-11 10:46:08,480 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:46:09,188 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:46:09,624 - tests.rls_test_helper - ERROR - Error creating test data: there is no unique or exclusion constraint matching the ON CONFLICT specification

2025-05-11 10:46:09,628 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:46:10,121 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:46:10,616 - tests.rls_test_helper - ERROR - Error creating test data: there is no unique or exclusion constraint matching the ON CONFLICT specification

2025-05-11 10:46:10,618 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:46:11,057 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:46:11,463 - tests.rls_test_helper - ERROR - Error creating test data: there is no unique or exclusion constraint matching the ON CONFLICT specification

2025-05-11 10:46:11,465 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:46:11,893 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:46:12,286 - tests.rls_test_helper - ERROR - Error creating test data: there is no unique or exclusion constraint matching the ON CONFLICT specification

2025-05-11 10:47:56,709 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:47:57,460 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:47:58,161 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:47:58,324 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:47:58,527 - tests.rls_test_helper - ERROR - Error creating test data: insert or update on table "team_members" violates foreign key constraint "team_members_user_id_fkey"
DETAIL:  Key (user_id)=(00000000-0000-0000-0000-000000000001) is not present in table "users".

2025-05-11 10:47:58,531 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:47:58,970 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:47:59,718 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:47:59,876 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:48:00,100 - tests.rls_test_helper - ERROR - Error creating test data: insert or update on table "team_members" violates foreign key constraint "team_members_user_id_fkey"
DETAIL:  Key (user_id)=(00000000-0000-0000-0000-000000000001) is not present in table "users".

2025-05-11 10:48:00,103 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:48:00,618 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:48:01,274 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:48:01,461 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:48:01,699 - tests.rls_test_helper - ERROR - Error creating test data: insert or update on table "team_members" violates foreign key constraint "team_members_user_id_fkey"
DETAIL:  Key (user_id)=(00000000-0000-0000-0000-000000000001) is not present in table "users".

2025-05-11 10:48:01,702 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:48:02,159 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:48:02,786 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:48:02,910 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:48:03,120 - tests.rls_test_helper - ERROR - Error creating test data: insert or update on table "team_members" violates foreign key constraint "team_members_user_id_fkey"
DETAIL:  Key (user_id)=(00000000-0000-0000-0000-000000000001) is not present in table "users".

2025-05-11 10:48:03,123 - __main__ - INFO - Report generated: reports/rls_verification_report_20250511_104803.md
2025-05-11 10:49:04,931 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:49:05,507 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:49:05,715 - tests.rls_test_helper - ERROR - Error creating test data: relation "users" does not exist
LINE 1: SELECT id FROM users WHERE id = '00000000-0000-0000-0000-000...
                       ^

2025-05-11 10:49:05,718 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:49:06,299 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:49:06,509 - tests.rls_test_helper - ERROR - Error creating test data: relation "users" does not exist
LINE 1: SELECT id FROM users WHERE id = '00000000-0000-0000-0000-000...
                       ^

2025-05-11 10:49:06,511 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:49:06,920 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:49:07,138 - tests.rls_test_helper - ERROR - Error creating test data: relation "users" does not exist
LINE 1: SELECT id FROM users WHERE id = '00000000-0000-0000-0000-000...
                       ^

2025-05-11 10:49:07,141 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:49:07,572 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:49:07,816 - tests.rls_test_helper - ERROR - Error creating test data: relation "users" does not exist
LINE 1: SELECT id FROM users WHERE id = '00000000-0000-0000-0000-000...
                       ^

2025-05-11 10:49:07,819 - __main__ - INFO - Report generated: reports/rls_verification_report_20250511_104907.md
2025-05-11 10:50:20,373 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:50:20,917 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:50:21,209 - tests.rls_test_helper - ERROR - Error creating test data: cannot insert a non-DEFAULT value into column "confirmed_at"
DETAIL:  Column "confirmed_at" is a generated column.

2025-05-11 10:50:21,212 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:50:21,619 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:50:21,901 - tests.rls_test_helper - ERROR - Error creating test data: cannot insert a non-DEFAULT value into column "confirmed_at"
DETAIL:  Column "confirmed_at" is a generated column.

2025-05-11 10:50:21,904 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:50:22,294 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:50:22,558 - tests.rls_test_helper - ERROR - Error creating test data: cannot insert a non-DEFAULT value into column "confirmed_at"
DETAIL:  Column "confirmed_at" is a generated column.

2025-05-11 10:50:22,561 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:50:22,963 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:50:23,240 - tests.rls_test_helper - ERROR - Error creating test data: cannot insert a non-DEFAULT value into column "confirmed_at"
DETAIL:  Column "confirmed_at" is a generated column.

2025-05-11 10:50:23,243 - __main__ - INFO - Report generated: reports/rls_verification_report_20250511_105023.md
2025-05-11 10:51:02,428 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:51:02,838 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:51:03,062 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 10:51:03,200 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 10:51:03,717 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:51:03,838 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:03,968 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:04,180 - tests.rls_test_helper - ERROR - Error creating test data: new row for relation "team_members" violates check constraint "team_members_role_check"
DETAIL:  Failing row contains (cd9f31b0-8911-4243-a793-04273d916a4f, 33333333-3333-3333-3333-333333333333, 00000000-0000-0000-0000-000000000002, member, 2025-05-11 16:51:02.867643+00, 2025-05-11 16:51:02.867643+00, null).

2025-05-11 10:51:04,182 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:51:04,577 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:51:04,781 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 10:51:04,917 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 10:51:05,758 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:51:05,920 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:06,059 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:06,248 - tests.rls_test_helper - ERROR - Error creating test data: new row for relation "team_members" violates check constraint "team_members_role_check"
DETAIL:  Failing row contains (ecc5801f-0f1f-4878-a009-cd8e8df9954f, 33333333-3333-3333-3333-333333333333, 00000000-0000-0000-0000-000000000002, member, 2025-05-11 16:51:04.605839+00, 2025-05-11 16:51:04.605839+00, null).

2025-05-11 10:51:06,251 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:51:06,719 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:51:07,006 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 10:51:07,155 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 10:51:07,808 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:51:07,999 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:08,219 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:08,416 - tests.rls_test_helper - ERROR - Error creating test data: new row for relation "team_members" violates check constraint "team_members_role_check"
DETAIL:  Failing row contains (e9b91b8e-9376-4ada-bdd1-46ec4c3bf605, 33333333-3333-3333-3333-333333333333, 00000000-0000-0000-0000-000000000002, member, 2025-05-11 16:51:06.748997+00, 2025-05-11 16:51:06.748997+00, null).

2025-05-11 10:51:08,419 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:51:08,852 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:51:09,116 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 10:51:09,258 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 10:51:09,869 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:51:10,004 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:10,156 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:10,360 - tests.rls_test_helper - ERROR - Error creating test data: new row for relation "team_members" violates check constraint "team_members_role_check"
DETAIL:  Failing row contains (25eb2eb3-726f-4c00-8cfa-4bca839d524a, 33333333-3333-3333-3333-333333333333, 00000000-0000-0000-0000-000000000002, member, 2025-05-11 16:51:08.907413+00, 2025-05-11 16:51:08.907413+00, null).

2025-05-11 10:51:10,364 - __main__ - INFO - Report generated: reports/rls_verification_report_20250511_105110.md
2025-05-11 10:51:31,153 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:51:31,797 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:51:32,041 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 10:51:32,173 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 10:51:32,786 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:51:32,927 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:33,099 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:33,303 - tests.rls_test_helper - INFO - Linked user2 00000000-0000-0000-0000-000000000002 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:33,588 - tests.rls_test_helper - ERROR - Error creating test data: column "project_id" does not exist
LINE 3:                 WHERE project_id IN ('11111111-1111-1111-111...
                              ^

2025-05-11 10:51:33,593 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:51:34,004 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:51:34,202 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 10:51:34,335 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 10:51:35,024 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:51:35,160 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:35,285 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:35,447 - tests.rls_test_helper - INFO - Linked user2 00000000-0000-0000-0000-000000000002 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:35,633 - tests.rls_test_helper - ERROR - Error creating test data: column "project_id" does not exist
LINE 3:                 WHERE project_id IN ('11111111-1111-1111-111...
                              ^

2025-05-11 10:51:35,636 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:51:36,030 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:51:36,250 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 10:51:36,385 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 10:51:36,957 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:51:37,139 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:37,297 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:37,431 - tests.rls_test_helper - INFO - Linked user2 00000000-0000-0000-0000-000000000002 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:37,689 - tests.rls_test_helper - ERROR - Error creating test data: column "project_id" does not exist
LINE 3:                 WHERE project_id IN ('11111111-1111-1111-111...
                              ^

2025-05-11 10:51:37,691 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 10:51:38,205 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 10:51:38,461 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 10:51:38,604 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 10:51:39,244 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 10:51:39,427 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:39,563 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:39,695 - tests.rls_test_helper - INFO - Linked user2 00000000-0000-0000-0000-000000000002 to team 33333333-3333-3333-3333-333333333333
2025-05-11 10:51:39,931 - tests.rls_test_helper - ERROR - Error creating test data: column "project_id" does not exist
LINE 3:                 WHERE project_id IN ('11111111-1111-1111-111...
                              ^

2025-05-11 10:51:39,934 - __main__ - INFO - Report generated: reports/rls_verification_report_20250511_105139.md
2025-05-11 11:03:10,072 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 11:03:10,689 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 11:03:10,919 - tests.rls_test_helper - INFO - Created test user1 in auth.users table with ID 00000000-0000-0000-0000-000000000001
2025-05-11 11:03:11,069 - tests.rls_test_helper - INFO - Created test user2 in auth.users table with ID 00000000-0000-0000-0000-000000000002
2025-05-11 11:03:11,706 - tests.rls_test_helper - INFO - Created 2 user profiles
2025-05-11 11:03:11,863 - tests.rls_test_helper - INFO - Created test team with ID 33333333-3333-3333-3333-333333333333
2025-05-11 11:03:12,187 - tests.rls_test_helper - INFO - Linked user1 00000000-0000-0000-0000-000000000001 to team 33333333-3333-3333-3333-333333333333
2025-05-11 11:03:12,319 - tests.rls_test_helper - INFO - Linked user2 00000000-0000-0000-0000-000000000002 to team 33333333-3333-3333-3333-333333333333
2025-05-11 11:03:12,647 - tests.rls_test_helper - INFO - Created test molecule 1 with ID ec8054d3-7813-4f22-9879-976e1e6b12ab
2025-05-11 11:03:12,720 - tests.rls_test_helper - INFO - Created test molecule 2 with ID f68cf7d8-16de-42e1-a07c-1e68f850157d
2025-05-11 11:03:12,821 - tests.rls_test_helper - INFO - Created experiment 1 linking molecule 1 to project 1
2025-05-11 11:03:12,898 - tests.rls_test_helper - INFO - Created experiment 2 linking molecule 2 to project 2
2025-05-11 11:03:13,062 - tests.rls_test_helper - INFO - Created test mixture 1 with ID 85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f
2025-05-11 11:03:13,130 - tests.rls_test_helper - INFO - Created test mixture 2 with ID 40cbb856-3b48-4b9e-841f-52fa27ac9ff1
2025-05-11 11:03:13,195 - tests.rls_test_helper - INFO - Test data created successfully
2025-05-11 11:03:13,401 - tests.test_rls_policies_verification - INFO - Test data created: {"users": ["00000000-0000-0000-0000-000000000001", "00000000-0000-0000-0000-000000000002"], "projects": ["11111111-1111-1111-1111-111111111111", "22222222-2222-2222-2222-222222222222"], "team": "33333333-3333-3333-3333-333333333333", "molecules": ["ec8054d3-7813-4f22-9879-976e1e6b12ab", "f68cf7d8-16de-42e1-a07c-1e68f850157d"], "experiments": ["de15deef-3c99-4e68-8710-5a67f70a76b1", "408fa489-6ea5-4f64-9167-fb04425ad66c"], "mixtures": ["85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f", "40cbb856-3b48-4b9e-841f-52fa27ac9ff1"]}
2025-05-11 11:03:14,008 - tests.rls_test_helper - INFO - Testing function is_project_member for user 00000000-0000-0000-0000-000000000001
2025-05-11 11:03:14,008 - tests.rls_test_helper - INFO - Executing function: SELECT is_project_member('22222222-2222-2222-2222-222222222222')
2025-05-11 11:03:14,102 - tests.rls_test_helper - ERROR - Error testing security definer function: function is_project_member(unknown) does not exist
LINE 1: SELECT is_project_member('22222222-2222-2222-2222-2222222222...
               ^
HINT:  No function matches the given name and argument types. You might need to add explicit type casts.

2025-05-11 11:03:14,811 - tests.rls_test_helper - INFO - Testing function is_project_member for user 00000000-0000-0000-0000-000000000001
2025-05-11 11:03:14,811 - tests.rls_test_helper - INFO - Executing function: SELECT is_project_member('11111111-1111-1111-1111-111111111111')
2025-05-11 11:03:14,877 - tests.rls_test_helper - ERROR - Error testing security definer function: function is_project_member(unknown) does not exist
LINE 1: SELECT is_project_member('11111111-1111-1111-1111-1111111111...
               ^
HINT:  No function matches the given name and argument types. You might need to add explicit type casts.

2025-05-11 11:03:15,565 - tests.rls_test_helper - INFO - Testing function is_project_member for user 00000000-0000-0000-0000-000000000002
2025-05-11 11:03:15,565 - tests.rls_test_helper - INFO - Executing function: SELECT is_project_member('22222222-2222-2222-2222-222222222222')
2025-05-11 11:03:15,632 - tests.rls_test_helper - ERROR - Error testing security definer function: function is_project_member(unknown) does not exist
LINE 1: SELECT is_project_member('22222222-2222-2222-2222-2222222222...
               ^
HINT:  No function matches the given name and argument types. You might need to add explicit type casts.

2025-05-11 11:03:16,297 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'is_project_owner'...
2025-05-11 11:03:16,369 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:16,957 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'is_team_member'...
2025-05-11 11:03:17,027 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:17,673 - tests.rls_test_helper - INFO - Testing function is_team_member for user 00000000-0000-0000-0000-000000000001
2025-05-11 11:03:17,673 - tests.rls_test_helper - INFO - Executing function: SELECT is_team_member('33333333-3333-3333-3333-333333333333')
2025-05-11 11:03:17,744 - tests.rls_test_helper - INFO - Function returned: False
2025-05-11 11:03:18,302 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'mixture_in_user_project'...
2025-05-11 11:03:18,365 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:18,962 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'molecule_in_user_project'...
2025-05-11 11:03:19,040 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:19,588 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'user_projects'...
2025-05-11 11:03:19,668 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:20,558 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): SELECT user_projects()...
2025-05-11 11:03:20,627 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:20,629 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 11:03:21,084 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 11:03:21,870 - tests.rls_test_helper - INFO - Test data created successfully
2025-05-11 11:03:22,313 - tests.test_rls_policies_verification - INFO - Test data created: {"users": ["00000000-0000-0000-0000-000000000001", "00000000-0000-0000-0000-000000000002"], "projects": ["11111111-1111-1111-1111-111111111111", "22222222-2222-2222-2222-222222222222"], "team": "33333333-3333-3333-3333-333333333333", "molecules": ["ec8054d3-7813-4f22-9879-976e1e6b12ab", "f68cf7d8-16de-42e1-a07c-1e68f850157d"], "experiments": ["de15deef-3c99-4e68-8710-5a67f70a76b1", "408fa489-6ea5-4f64-9167-fb04425ad66c"], "mixtures": ["85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f", "40cbb856-3b48-4b9e-841f-52fa27ac9ff1"]}
2025-05-11 11:03:22,890 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('61ef852e-0495...
2025-05-11 11:03:22,978 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:23,567 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
           ...
2025-05-11 11:03:23,638 - tests.rls_test_helper - ERROR - Error executing SQL: insert or update on table "experiments" violates foreign key constraint "experiments_molecule_id_fkey"
DETAIL:  Key (molecule_id)=(61ef852e-0495-4a44-9aab-ad4cd7f47f83) is not present in table "molecules".

2025-05-11 11:03:24,220 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            DELETE FROM experiments
            WHERE molecule_id = '61ef852e-0495-4a44-9aab-ad4cd7...
2025-05-11 11:03:24,290 - tests.rls_test_helper - INFO - Query affected -1 rows
2025-05-11 11:03:25,290 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            DELETE FROM molecules
            WHERE id = '61ef852e-0495-4a44-9aab-ad4cd7f47f83'
   ...
2025-05-11 11:03:25,355 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:03:25,968 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('1018361d-b69d...
2025-05-11 11:03:26,054 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:26,678 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
           ...
2025-05-11 11:03:26,748 - tests.rls_test_helper - ERROR - Error executing SQL: insert or update on table "experiments" violates foreign key constraint "experiments_molecule_id_fkey"
DETAIL:  Key (molecule_id)=(1018361d-b69d-4885-87f6-f8f4170c19eb) is not present in table "molecules".

2025-05-11 11:03:27,399 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            DELETE FROM molecules
            WHERE id = '1018361d-b69d-4885-87f6-f8f4170c19eb'
   ...
2025-05-11 11:03:27,478 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:03:28,173 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('e3fb7982-4501...
2025-05-11 11:03:28,239 - tests.rls_test_helper - ERROR - Error executing SQL: new row violates row-level security policy for table "molecules"

2025-05-11 11:03:28,886 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('f9bd9c88-9b6d...
2025-05-11 11:03:28,966 - tests.rls_test_helper - ERROR - Error executing SQL: new row violates row-level security policy for table "molecules"

2025-05-11 11:03:29,704 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to molecules
2025-05-11 11:03:29,704 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM molecules WHERE id IN (SELECT molecule_id FROM experiments WHERE project_id = '11111111-1111-1111-1111-111111111111')
2025-05-11 11:03:29,805 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:03:30,756 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to molecules
2025-05-11 11:03:30,756 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM molecules WHERE id IN (SELECT molecule_id FROM experiments WHERE project_id = '22222222-2222-2222-2222-222222222222')
2025-05-11 11:03:30,852 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:03:31,921 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            SELECT m.id
            FROM molecules m
            JOIN experiments e ON m.id = e.mol...
2025-05-11 11:03:31,989 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:32,660 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            UPDATE molecules
            SET name = 'Updated Test Molecule'
            WHERE id = ...
2025-05-11 11:03:32,727 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:03:33,268 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            SELECT m.id
            FROM molecules m
            JOIN experiments e ON m.id = e.mol...
2025-05-11 11:03:33,338 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:34,083 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            UPDATE molecules
            SET name = 'Should Fail Update'
            WHERE id = 'f6...
2025-05-11 11:03:34,150 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:03:34,152 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 11:03:34,595 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 11:03:35,556 - tests.rls_test_helper - INFO - Test data created successfully
2025-05-11 11:03:36,225 - tests.test_rls_policies_verification - INFO - Test data created: {"users": ["00000000-0000-0000-0000-000000000001", "00000000-0000-0000-0000-000000000002"], "projects": ["11111111-1111-1111-1111-111111111111", "22222222-2222-2222-2222-222222222222"], "team": "33333333-3333-3333-3333-333333333333", "molecules": ["ec8054d3-7813-4f22-9879-976e1e6b12ab", "f68cf7d8-16de-42e1-a07c-1e68f850157d"], "experiments": ["de15deef-3c99-4e68-8710-5a67f70a76b1", "408fa489-6ea5-4f64-9167-fb04425ad66c"], "mixtures": ["85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f", "40cbb856-3b48-4b9e-841f-52fa27ac9ff1"]}
2025-05-11 11:03:36,792 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create a molecule for project1
            WITH mol1 AS (
                SEL...
2025-05-11 11:03:36,896 - tests.rls_test_helper - ERROR - Error executing SQL: column "project_id" does not exist
LINE 4:                 SELECT id FROM molecules WHERE project_id = ...
                                                       ^

2025-05-11 11:03:37,567 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            WITH mol2 AS (
                SELECT id FROM molecules WHERE project_id = '22222222-22...
2025-05-11 11:03:37,641 - tests.rls_test_helper - ERROR - Error executing SQL: column "project_id" does not exist
LINE 3:                 SELECT id FROM molecules WHERE project_id = ...
                                                       ^

2025-05-11 11:03:38,438 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create a mixture for project1
            WITH mix1 AS (
                SELE...
2025-05-11 11:03:38,516 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:39,250 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO mixture_components (mixture_id, molecule_id, concentration)
           ...
2025-05-11 11:03:39,337 - tests.rls_test_helper - ERROR - Error executing SQL: null value in column "concentration_unit" of relation "mixture_components" violates not-null constraint
DETAIL:  Failing row contains (1c449a68-b180-4d48-a594-591abdd139b1, 85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f, ec8054d3-7813-4f22-9879-976e1e6b12ab, 50.0, null, null, null, 2025-05-11 17:03:39.145507+00, 2025-05-11 17:03:39.145507+00, {}, null).

2025-05-11 11:03:39,920 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create a mixture for project2
            WITH mix2 AS (
                SELE...
2025-05-11 11:03:39,987 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:40,579 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO mixture_components (mixture_id, molecule_id, concentration)
           ...
2025-05-11 11:03:40,644 - tests.rls_test_helper - ERROR - Error executing SQL: null value in column "concentration_unit" of relation "mixture_components" violates not-null constraint
DETAIL:  Failing row contains (241ec126-0897-47cb-9f1e-a0f004063a98, 40cbb856-3b48-4b9e-841f-52fa27ac9ff1, f68cf7d8-16de-42e1-a07c-1e68f850157d, 75.0, null, null, null, 2025-05-11 17:03:40.474608+00, 2025-05-11 17:03:40.474608+00, {}, null).

2025-05-11 11:03:41,242 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create an experiment for project1
            WITH exp1 AS (
                ...
2025-05-11 11:03:41,304 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:41,868 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO experiment_properties (experiment_id, property_name, property_value, pr...
2025-05-11 11:03:41,933 - tests.rls_test_helper - ERROR - Error executing SQL: column "property_name" of relation "experiment_properties" does not exist
LINE 2: ...INSERT INTO experiment_properties (experiment_id, property_n...
                                                             ^

2025-05-11 11:03:42,520 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create an experiment for project2
            WITH exp2 AS (
                ...
2025-05-11 11:03:42,663 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:43,263 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO experiment_properties (experiment_id, property_name, property_value, pr...
2025-05-11 11:03:43,330 - tests.rls_test_helper - ERROR - Error executing SQL: column "property_name" of relation "experiment_properties" does not exist
LINE 2: ...INSERT INTO experiment_properties (experiment_id, property_n...
                                                             ^

2025-05-11 11:03:44,018 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to experiment_properties
2025-05-11 11:03:44,018 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM experiment_properties WHERE experiment_id = '408fa489-6ea5-4f64-9167-fb04425ad66c'
2025-05-11 11:03:44,144 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:03:44,896 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to experiment_properties
2025-05-11 11:03:44,896 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM experiment_properties WHERE experiment_id = 'de15deef-3c99-4e68-8710-5a67f70a76b1'
2025-05-11 11:03:45,010 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:03:45,795 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to mixture_components
2025-05-11 11:03:45,795 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM mixture_components WHERE mixture_id = '40cbb856-3b48-4b9e-841f-52fa27ac9ff1'
2025-05-11 11:03:45,864 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:03:46,510 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to mixture_components
2025-05-11 11:03:46,510 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM mixture_components WHERE mixture_id = '85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f'
2025-05-11 11:03:46,611 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:03:47,243 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to molecular_properties
2025-05-11 11:03:47,243 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM molecular_properties WHERE molecule_id = 'f68cf7d8-16de-42e1-a07c-1e68f850157d'
2025-05-11 11:03:47,313 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:03:47,934 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to molecular_properties
2025-05-11 11:03:47,934 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM molecular_properties WHERE molecule_id = 'ec8054d3-7813-4f22-9879-976e1e6b12ab'
2025-05-11 11:03:48,001 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:03:48,003 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 11:03:48,408 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 11:03:49,425 - tests.rls_test_helper - INFO - Test data created successfully
2025-05-11 11:03:50,084 - tests.test_rls_policies_verification - INFO - Test data created: {"users": ["00000000-0000-0000-0000-000000000001", "00000000-0000-0000-0000-000000000002"], "projects": ["11111111-1111-1111-1111-111111111111", "22222222-2222-2222-2222-222222222222"], "team": "33333333-3333-3333-3333-333333333333", "molecules": ["ec8054d3-7813-4f22-9879-976e1e6b12ab", "f68cf7d8-16de-42e1-a07c-1e68f850157d"], "experiments": ["de15deef-3c99-4e68-8710-5a67f70a76b1", "408fa489-6ea5-4f64-9167-fb04425ad66c"], "mixtures": ["85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f", "40cbb856-3b48-4b9e-841f-52fa27ac9ff1"]}
2025-05-11 11:03:50,664 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO molecules (id, project_id, name, smiles)
            VALUES ('4df10da3-75d1...
2025-05-11 11:03:50,734 - tests.rls_test_helper - ERROR - Error executing SQL: column "project_id" of relation "molecules" does not exist
LINE 2:             INSERT INTO molecules (id, project_id, name, smi...
                                               ^

2025-05-11 11:03:51,443 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO molecules (id, project_id, name, smiles)
            VALUES ('fad20f63-3239...
2025-05-11 11:03:51,506 - tests.rls_test_helper - ERROR - Error executing SQL: column "project_id" of relation "molecules" does not exist
LINE 2:             INSERT INTO molecules (id, project_id, name, smi...
                                               ^

2025-05-11 11:03:52,469 - tests.rls_test_helper - INFO - Testing service role access to molecules
2025-05-11 11:03:52,543 - tests.rls_test_helper - INFO - Query returned 10 results
2025-05-11 11:03:53,164 - tests.rls_test_helper - INFO - Testing service role access to mixtures
2025-05-11 11:03:53,238 - tests.rls_test_helper - INFO - Query returned 10 results
2025-05-11 11:03:53,853 - tests.rls_test_helper - INFO - Testing service role access to experiments
2025-05-11 11:03:53,939 - tests.rls_test_helper - INFO - Query returned 5 results
2025-05-11 11:03:54,646 - tests.rls_test_helper - INFO - Testing service role access to projects
2025-05-11 11:03:54,716 - tests.rls_test_helper - INFO - Query returned 3 results
2025-05-11 11:03:55,298 - tests.rls_test_helper - INFO - Testing service role access to user_profile
2025-05-11 11:03:55,386 - tests.rls_test_helper - INFO - Query returned 3 results
2025-05-11 11:03:56,060 - tests.rls_test_helper - INFO - Testing service role access to molecular_properties
2025-05-11 11:03:56,155 - tests.rls_test_helper - INFO - Query returned 10 results
2025-05-11 11:03:56,704 - tests.rls_test_helper - INFO - Testing service role access to mixture_components
2025-05-11 11:03:56,794 - tests.rls_test_helper - INFO - Query returned 3 results
2025-05-11 11:03:57,690 - tests.rls_test_helper - INFO - Testing service role access to experiment_properties
2025-05-11 11:03:57,757 - tests.rls_test_helper - INFO - Query returned 3 results
2025-05-11 11:03:58,339 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT id FROM molecules LIMIT 1...
2025-05-11 11:03:58,408 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:58,931 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                UPDATE molecules
                SET name = 'Service Role Updated Test Molecule'
  ...
2025-05-11 11:03:59,012 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:03:59,014 - __main__ - INFO - Report generated: reports/rls_verification_report_20250511_110359.md
2025-05-11 11:06:28,958 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 11:06:29,585 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 11:06:30,467 - tests.rls_test_helper - INFO - Test data created successfully
2025-05-11 11:06:31,141 - tests.test_rls_policies_verification - INFO - Test data created: {"users": ["00000000-0000-0000-0000-000000000001", "00000000-0000-0000-0000-000000000002"], "projects": ["11111111-1111-1111-1111-111111111111", "22222222-2222-2222-2222-222222222222"], "team": "33333333-3333-3333-3333-333333333333", "molecules": ["ec8054d3-7813-4f22-9879-976e1e6b12ab", "f68cf7d8-16de-42e1-a07c-1e68f850157d"], "experiments": ["de15deef-3c99-4e68-8710-5a67f70a76b1", "408fa489-6ea5-4f64-9167-fb04425ad66c"], "mixtures": ["85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f", "40cbb856-3b48-4b9e-841f-52fa27ac9ff1"]}
2025-05-11 11:06:31,861 - tests.rls_test_helper - INFO - Testing function is_project_member for user 00000000-0000-0000-0000-000000000001
2025-05-11 11:06:31,861 - tests.rls_test_helper - INFO - Executing function: SELECT is_project_member('22222222-2222-2222-2222-222222222222')
2025-05-11 11:06:32,009 - tests.rls_test_helper - ERROR - Error testing security definer function: function is_project_member(unknown) does not exist
LINE 1: SELECT is_project_member('22222222-2222-2222-2222-2222222222...
               ^
HINT:  No function matches the given name and argument types. You might need to add explicit type casts.

2025-05-11 11:06:32,713 - tests.rls_test_helper - INFO - Testing function is_project_member for user 00000000-0000-0000-0000-000000000001
2025-05-11 11:06:32,714 - tests.rls_test_helper - INFO - Executing function: SELECT is_project_member('11111111-1111-1111-1111-111111111111')
2025-05-11 11:06:32,779 - tests.rls_test_helper - ERROR - Error testing security definer function: function is_project_member(unknown) does not exist
LINE 1: SELECT is_project_member('11111111-1111-1111-1111-1111111111...
               ^
HINT:  No function matches the given name and argument types. You might need to add explicit type casts.

2025-05-11 11:06:33,540 - tests.rls_test_helper - INFO - Testing function is_project_member for user 00000000-0000-0000-0000-000000000002
2025-05-11 11:06:33,541 - tests.rls_test_helper - INFO - Executing function: SELECT is_project_member('22222222-2222-2222-2222-222222222222')
2025-05-11 11:06:33,616 - tests.rls_test_helper - ERROR - Error testing security definer function: function is_project_member(unknown) does not exist
LINE 1: SELECT is_project_member('22222222-2222-2222-2222-2222222222...
               ^
HINT:  No function matches the given name and argument types. You might need to add explicit type casts.

2025-05-11 11:06:34,260 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'is_project_owner'...
2025-05-11 11:06:34,372 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:35,176 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'is_team_member'...
2025-05-11 11:06:35,244 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:35,852 - tests.rls_test_helper - INFO - Testing function is_team_member for user 00000000-0000-0000-0000-000000000001
2025-05-11 11:06:35,852 - tests.rls_test_helper - INFO - Executing function: SELECT is_team_member('33333333-3333-3333-3333-333333333333')
2025-05-11 11:06:35,920 - tests.rls_test_helper - INFO - Function returned: False
2025-05-11 11:06:36,594 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'mixture_in_user_project'...
2025-05-11 11:06:36,712 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:37,324 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'molecule_in_user_project'...
2025-05-11 11:06:37,395 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:37,995 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM pg_proc WHERE proname = 'user_projects'...
2025-05-11 11:06:38,060 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:38,698 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): SELECT user_projects()...
2025-05-11 11:06:38,769 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:38,770 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 11:06:39,169 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 11:06:40,114 - tests.rls_test_helper - INFO - Test data created successfully
2025-05-11 11:06:40,722 - tests.test_rls_policies_verification - INFO - Test data created: {"users": ["00000000-0000-0000-0000-000000000001", "00000000-0000-0000-0000-000000000002"], "projects": ["11111111-1111-1111-1111-111111111111", "22222222-2222-2222-2222-222222222222"], "team": "33333333-3333-3333-3333-333333333333", "molecules": ["ec8054d3-7813-4f22-9879-976e1e6b12ab", "f68cf7d8-16de-42e1-a07c-1e68f850157d"], "experiments": ["de15deef-3c99-4e68-8710-5a67f70a76b1", "408fa489-6ea5-4f64-9167-fb04425ad66c"], "mixtures": ["85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f", "40cbb856-3b48-4b9e-841f-52fa27ac9ff1"]}
2025-05-11 11:06:41,385 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('688ca30d-8ff8...
2025-05-11 11:06:41,512 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:42,122 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
           ...
2025-05-11 11:06:42,194 - tests.rls_test_helper - ERROR - Error executing SQL: insert or update on table "experiments" violates foreign key constraint "experiments_molecule_id_fkey"
DETAIL:  Key (molecule_id)=(688ca30d-8ff8-48a9-86e4-37fb7e010b46) is not present in table "molecules".

2025-05-11 11:06:43,240 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            DELETE FROM experiments
            WHERE molecule_id = '688ca30d-8ff8-48a9-86e4-37fb7e...
2025-05-11 11:06:43,323 - tests.rls_test_helper - INFO - Query affected -1 rows
2025-05-11 11:06:44,353 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            DELETE FROM molecules
            WHERE id = '688ca30d-8ff8-48a9-86e4-37fb7e010b46'
   ...
2025-05-11 11:06:45,051 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:06:45,623 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('1458a1ea-a336...
2025-05-11 11:06:45,710 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:46,415 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO experiments (id, molecule_id, project_id, name, created_by)
           ...
2025-05-11 11:06:46,488 - tests.rls_test_helper - ERROR - Error executing SQL: insert or update on table "experiments" violates foreign key constraint "experiments_molecule_id_fkey"
DETAIL:  Key (molecule_id)=(1458a1ea-a336-4dd1-8d57-b9ac5da7f136) is not present in table "molecules".

2025-05-11 11:06:47,396 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            DELETE FROM molecules
            WHERE id = '1458a1ea-a336-4dd1-8d57-b9ac5da7f136'
   ...
2025-05-11 11:06:47,566 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:06:48,199 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('641a093e-203c...
2025-05-11 11:06:48,270 - tests.rls_test_helper - ERROR - Error executing SQL: new row violates row-level security policy for table "molecules"

2025-05-11 11:06:48,987 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('19ebc790-bfcf...
2025-05-11 11:06:49,124 - tests.rls_test_helper - ERROR - Error executing SQL: new row violates row-level security policy for table "molecules"

2025-05-11 11:06:49,885 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to molecules
2025-05-11 11:06:49,885 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM molecules WHERE id IN (SELECT molecule_id FROM experiments WHERE project_id = '11111111-1111-1111-1111-111111111111')
2025-05-11 11:06:49,952 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:06:50,758 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to molecules
2025-05-11 11:06:50,758 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM molecules WHERE id IN (SELECT molecule_id FROM experiments WHERE project_id = '22222222-2222-2222-2222-222222222222')
2025-05-11 11:06:50,825 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:06:51,531 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            SELECT m.id
            FROM molecules m
            JOIN experiments e ON m.id = e.mol...
2025-05-11 11:06:51,680 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:52,371 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            UPDATE molecules
            SET name = 'Updated Test Molecule'
            WHERE id = ...
2025-05-11 11:06:52,435 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:06:53,182 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            SELECT m.id
            FROM molecules m
            JOIN experiments e ON m.id = e.mol...
2025-05-11 11:06:53,254 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:53,913 - tests.rls_test_helper - INFO - Executing SQL as authenticated (user_id: 00000000-0000-0000-0000-000000000001): 
            UPDATE molecules
            SET name = 'Should Fail Update'
            WHERE id = 'f6...
2025-05-11 11:06:53,996 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:06:53,997 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 11:06:54,431 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 11:06:55,530 - tests.rls_test_helper - INFO - Test data created successfully
2025-05-11 11:06:56,157 - tests.test_rls_policies_verification - INFO - Test data created: {"users": ["00000000-0000-0000-0000-000000000001", "00000000-0000-0000-0000-000000000002"], "projects": ["11111111-1111-1111-1111-111111111111", "22222222-2222-2222-2222-222222222222"], "team": "33333333-3333-3333-3333-333333333333", "molecules": ["ec8054d3-7813-4f22-9879-976e1e6b12ab", "f68cf7d8-16de-42e1-a07c-1e68f850157d"], "experiments": ["de15deef-3c99-4e68-8710-5a67f70a76b1", "408fa489-6ea5-4f64-9167-fb04425ad66c"], "mixtures": ["85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f", "40cbb856-3b48-4b9e-841f-52fa27ac9ff1"]}
2025-05-11 11:06:56,767 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM molecular_properties WHERE molecule_id = 'ec8054d3-7813-4f22-9879-976e1e6b12ab'...
2025-05-11 11:06:56,831 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:57,491 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                    INSERT INTO molecular_properties (id, molecule_id, property_name, property_valu...
2025-05-11 11:06:57,568 - tests.rls_test_helper - ERROR - Error executing SQL: null value in column "property_type_id" of relation "molecular_properties" violates not-null constraint
DETAIL:  Failing row contains (bd642416-4904-4d41-8af5-19e4edd1c0f2, ec8054d3-7813-4f22-9879-976e1e6b12ab, null, null, null, null, C, 00000000-0000-0000-0000-000000000001, null, null, null, 2025-05-11 17:06:57.362806+00, 2025-05-11 17:06:57.362806+00, null, boiling_point, 100.0, null, null).

2025-05-11 11:06:58,214 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM molecular_properties WHERE molecule_id = 'f68cf7d8-16de-42e1-a07c-1e68f850157d'...
2025-05-11 11:06:58,302 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:06:59,158 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                    INSERT INTO molecular_properties (id, molecule_id, property_name, property_valu...
2025-05-11 11:06:59,229 - tests.rls_test_helper - ERROR - Error executing SQL: null value in column "property_type_id" of relation "molecular_properties" violates not-null constraint
DETAIL:  Failing row contains (4e1bc3e1-7160-45a1-b438-27625c4db294, f68cf7d8-16de-42e1-a07c-1e68f850157d, null, null, null, null, C, 00000000-0000-0000-0000-000000000002, null, null, null, 2025-05-11 17:06:58.866876+00, 2025-05-11 17:06:58.866876+00, null, boiling_point, 120.0, null, null).

2025-05-11 11:06:59,752 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM mixture_components WHERE mixture_id = '85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f'...
2025-05-11 11:06:59,817 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:00,947 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                    INSERT INTO mixture_components (id, mixture_id, molecule_id, concentration, con...
2025-05-11 11:07:01,080 - tests.rls_test_helper - INFO - Query affected -1 rows
2025-05-11 11:07:01,800 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM mixture_components WHERE mixture_id = '40cbb856-3b48-4b9e-841f-52fa27ac9ff1'...
2025-05-11 11:07:01,887 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:02,438 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                    INSERT INTO mixture_components (id, mixture_id, molecule_id, concentration, con...
2025-05-11 11:07:02,509 - tests.rls_test_helper - INFO - Query affected -1 rows
2025-05-11 11:07:03,175 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM experiment_properties WHERE experiment_id = 'de15deef-3c99-4e68-8710-5a67f70a76...
2025-05-11 11:07:03,288 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:04,149 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                    INSERT INTO experiment_properties (id, experiment_id, property_type_id, numeric...
2025-05-11 11:07:04,306 - tests.rls_test_helper - ERROR - Error executing SQL: insert or update on table "experiment_properties" violates foreign key constraint "experiment_properties_property_type_id_fkey"
DETAIL:  Key (property_type_id)=(8c234177-8362-473e-9807-5ad81ed687fa) is not present in table "property_types".

2025-05-11 11:07:04,972 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT COUNT(*) FROM experiment_properties WHERE experiment_id = '408fa489-6ea5-4f64-9167-fb04425ad6...
2025-05-11 11:07:05,044 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:05,683 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                    INSERT INTO experiment_properties (id, experiment_id, property_type_id, numeric...
2025-05-11 11:07:05,760 - tests.rls_test_helper - ERROR - Error executing SQL: insert or update on table "experiment_properties" violates foreign key constraint "experiment_properties_property_type_id_fkey"
DETAIL:  Key (property_type_id)=(e13e0eee-68b5-483d-ab0a-a1802e0137ca) is not present in table "property_types".

2025-05-11 11:07:06,343 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create a molecule for project1
            WITH mol1 AS (
                SEL...
2025-05-11 11:07:06,419 - tests.rls_test_helper - ERROR - Error executing SQL: column "project_id" does not exist
LINE 4:                 SELECT id FROM molecules WHERE project_id = ...
                                                       ^

2025-05-11 11:07:07,231 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            WITH mol2 AS (
                SELECT id FROM molecules WHERE project_id = '22222222-22...
2025-05-11 11:07:07,316 - tests.rls_test_helper - ERROR - Error executing SQL: column "project_id" does not exist
LINE 3:                 SELECT id FROM molecules WHERE project_id = ...
                                                       ^

2025-05-11 11:07:07,843 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create a mixture for project1
            WITH mix1 AS (
                SELE...
2025-05-11 11:07:07,913 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:08,486 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO mixture_components (mixture_id, molecule_id, concentration)
           ...
2025-05-11 11:07:08,562 - tests.rls_test_helper - ERROR - Error executing SQL: null value in column "concentration_unit" of relation "mixture_components" violates not-null constraint
DETAIL:  Failing row contains (42b5061e-f193-41f7-9c13-b21be72d4951, 85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f, ec8054d3-7813-4f22-9879-976e1e6b12ab, 50.0, null, null, null, 2025-05-11 17:07:08.387102+00, 2025-05-11 17:07:08.387102+00, {}, null).

2025-05-11 11:07:09,230 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create a mixture for project2
            WITH mix2 AS (
                SELE...
2025-05-11 11:07:09,384 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:10,021 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO mixture_components (mixture_id, molecule_id, concentration)
           ...
2025-05-11 11:07:10,086 - tests.rls_test_helper - ERROR - Error executing SQL: null value in column "concentration_unit" of relation "mixture_components" violates not-null constraint
DETAIL:  Failing row contains (ff66db9f-113c-4f09-ac5c-655a08cf9211, 40cbb856-3b48-4b9e-841f-52fa27ac9ff1, f68cf7d8-16de-42e1-a07c-1e68f850157d, 75.0, null, null, null, 2025-05-11 17:07:09.9236+00, 2025-05-11 17:07:09.9236+00, {}, null).

2025-05-11 11:07:10,765 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create an experiment for project1
            WITH exp1 AS (
                ...
2025-05-11 11:07:10,832 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:11,534 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO experiment_properties (experiment_id, property_name, property_value, pr...
2025-05-11 11:07:11,602 - tests.rls_test_helper - ERROR - Error executing SQL: column "property_name" of relation "experiment_properties" does not exist
LINE 2: ...INSERT INTO experiment_properties (experiment_id, property_n...
                                                             ^

2025-05-11 11:07:12,416 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            -- Get or create an experiment for project2
            WITH exp2 AS (
                ...
2025-05-11 11:07:12,483 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:13,324 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                INSERT INTO experiment_properties (experiment_id, property_name, property_value, pr...
2025-05-11 11:07:13,399 - tests.rls_test_helper - ERROR - Error executing SQL: column "property_name" of relation "experiment_properties" does not exist
LINE 2: ...INSERT INTO experiment_properties (experiment_id, property_n...
                                                             ^

2025-05-11 11:07:14,090 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to experiment_properties
2025-05-11 11:07:14,091 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM experiment_properties WHERE experiment_id = '408fa489-6ea5-4f64-9167-fb04425ad66c'
2025-05-11 11:07:14,168 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:07:15,199 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to experiment_properties
2025-05-11 11:07:15,200 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM experiment_properties WHERE experiment_id = 'de15deef-3c99-4e68-8710-5a67f70a76b1'
2025-05-11 11:07:15,344 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:07:16,884 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to mixture_components
2025-05-11 11:07:16,885 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM mixture_components WHERE mixture_id = '40cbb856-3b48-4b9e-841f-52fa27ac9ff1'
2025-05-11 11:07:16,960 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:07:17,652 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to mixture_components
2025-05-11 11:07:17,652 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM mixture_components WHERE mixture_id = '85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f'
2025-05-11 11:07:17,718 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:07:18,382 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to molecular_properties
2025-05-11 11:07:18,382 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM molecular_properties WHERE molecule_id = 'f68cf7d8-16de-42e1-a07c-1e68f850157d'
2025-05-11 11:07:18,449 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:07:19,146 - tests.rls_test_helper - INFO - Testing user 00000000-0000-0000-0000-000000000001 access to molecular_properties
2025-05-11 11:07:19,146 - tests.rls_test_helper - INFO - Executing query: SELECT * FROM molecular_properties WHERE molecule_id = 'ec8054d3-7813-4f22-9879-976e1e6b12ab'
2025-05-11 11:07:19,228 - tests.rls_test_helper - ERROR - Error testing user access: infinite recursion detected in policy for relation "user_profile"

2025-05-11 11:07:19,230 - tests.rls_test_helper - INFO - RLSTestHelper initialized with connection parameters
2025-05-11 11:07:19,655 - tests.rls_test_helper - INFO - Creating test data for RLS testing
2025-05-11 11:07:20,511 - tests.rls_test_helper - INFO - Test data created successfully
2025-05-11 11:07:21,450 - tests.test_rls_policies_verification - INFO - Test data created: {"users": ["00000000-0000-0000-0000-000000000001", "00000000-0000-0000-0000-000000000002"], "projects": ["11111111-1111-1111-1111-111111111111", "22222222-2222-2222-2222-222222222222"], "team": "33333333-3333-3333-3333-333333333333", "molecules": ["ec8054d3-7813-4f22-9879-976e1e6b12ab", "f68cf7d8-16de-42e1-a07c-1e68f850157d"], "experiments": ["de15deef-3c99-4e68-8710-5a67f70a76b1", "408fa489-6ea5-4f64-9167-fb04425ad66c"], "mixtures": ["85f6a0dd-8ae2-4876-8c37-ed58fd0ca36f", "40cbb856-3b48-4b9e-841f-52fa27ac9ff1"]}
2025-05-11 11:07:22,059 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('6e156597-a55c...
2025-05-11 11:07:22,147 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:22,777 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            DELETE FROM molecules
            WHERE id = '6e156597-a55c-4bd7-8a9d-6270eee2acd1'
   ...
2025-05-11 11:07:22,848 - tests.rls_test_helper - INFO - Query returned 0 results
2025-05-11 11:07:23,475 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO molecules (id, name, smiles, created_by)
            VALUES ('b816d5ae-3f4a...
2025-05-11 11:07:23,614 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:24,248 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
            INSERT INTO mixtures (id, project_id, name, created_by)
            VALUES ('c05c8190-9...
2025-05-11 11:07:24,329 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:24,930 - tests.rls_test_helper - INFO - Testing service role access to molecules
2025-05-11 11:07:25,003 - tests.rls_test_helper - INFO - Query returned 10 results
2025-05-11 11:07:25,638 - tests.rls_test_helper - INFO - Testing service role access to mixtures
2025-05-11 11:07:25,767 - tests.rls_test_helper - INFO - Query returned 10 results
2025-05-11 11:07:26,445 - tests.rls_test_helper - INFO - Testing service role access to experiments
2025-05-11 11:07:26,522 - tests.rls_test_helper - INFO - Query returned 5 results
2025-05-11 11:07:27,153 - tests.rls_test_helper - INFO - Testing service role access to projects
2025-05-11 11:07:27,219 - tests.rls_test_helper - INFO - Query returned 3 results
2025-05-11 11:07:27,818 - tests.rls_test_helper - INFO - Testing service role access to user_profile
2025-05-11 11:07:27,910 - tests.rls_test_helper - INFO - Query returned 3 results
2025-05-11 11:07:28,514 - tests.rls_test_helper - INFO - Testing service role access to molecular_properties
2025-05-11 11:07:28,583 - tests.rls_test_helper - INFO - Query returned 10 results
2025-05-11 11:07:29,134 - tests.rls_test_helper - INFO - Testing service role access to mixture_components
2025-05-11 11:07:29,206 - tests.rls_test_helper - INFO - Query returned 3 results
2025-05-11 11:07:29,962 - tests.rls_test_helper - INFO - Testing service role access to experiment_properties
2025-05-11 11:07:30,049 - tests.rls_test_helper - INFO - Query returned 3 results
2025-05-11 11:07:30,670 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): SELECT id FROM molecules LIMIT 1...
2025-05-11 11:07:30,739 - tests.rls_test_helper - INFO - Query returned 1 results
2025-05-11 11:07:31,294 - tests.rls_test_helper - INFO - Executing SQL as service_role (user_id: None): 
                UPDATE molecules
                SET name = 'Service Role Updated Test Molecule'
  ...
2025-05-11 11:07:31,369 - tests.rls_test_helper - INFO - Query returned 1 results
```
