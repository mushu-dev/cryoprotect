# RLS Verification Report

Generated on: 2025-05-11 10:51:10

## Test Summary

- Tests run: 0
- Failures: 0
- Errors: 4
- Skipped: 0
- Success rate: N/A (no tests were run)

## Errors

### setUpClass (tests.test_rls_policies_verification.RLSSecurityDefinerFunctionTests)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 46, in setUpClass
    cls.test_data = cls.helper.create_test_data()
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 215, in create_test_data
    raise e
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 152, in create_test_data
    cursor.execute(f"""
    ~~~~~~~~~~~~~~^^^^^
        INSERT INTO team_members (team_id, user_id, role)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        VALUES ('{team_id}', '{user2_id}', 'member')
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    """)
    ^^^^
psycopg2.errors.CheckViolation: new row for relation "team_members" violates check constraint "team_members_role_check"
DETAIL:  Failing row contains (cd9f31b0-8911-4243-a793-04273d916a4f, 33333333-3333-3333-3333-333333333333, 00000000-0000-0000-0000-000000000002, member, 2025-05-11 16:51:02.867643+00, 2025-05-11 16:51:02.867643+00, null).

```

### setUpClass (tests.test_rls_policies_verification.RLSTableAccessPolicyTests)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 231, in setUpClass
    cls.test_data = cls.helper.create_test_data()
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 215, in create_test_data
    raise e
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 152, in create_test_data
    cursor.execute(f"""
    ~~~~~~~~~~~~~~^^^^^
        INSERT INTO team_members (team_id, user_id, role)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        VALUES ('{team_id}', '{user2_id}', 'member')
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    """)
    ^^^^
psycopg2.errors.CheckViolation: new row for relation "team_members" violates check constraint "team_members_role_check"
DETAIL:  Failing row contains (ecc5801f-0f1f-4878-a009-cd8e8df9954f, 33333333-3333-3333-3333-333333333333, 00000000-0000-0000-0000-000000000002, member, 2025-05-11 16:51:04.605839+00, 2025-05-11 16:51:04.605839+00, null).

```

### setUpClass (tests.test_rls_policies_verification.RLSRelationshipPolicyTests)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 380, in setUpClass
    cls.test_data = cls.helper.create_test_data()
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 215, in create_test_data
    raise e
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 152, in create_test_data
    cursor.execute(f"""
    ~~~~~~~~~~~~~~^^^^^
        INSERT INTO team_members (team_id, user_id, role)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        VALUES ('{team_id}', '{user2_id}', 'member')
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    """)
    ^^^^
psycopg2.errors.CheckViolation: new row for relation "team_members" violates check constraint "team_members_role_check"
DETAIL:  Failing row contains (e9b91b8e-9376-4ada-bdd1-46ec4c3bf605, 33333333-3333-3333-3333-333333333333, 00000000-0000-0000-0000-000000000002, member, 2025-05-11 16:51:06.748997+00, 2025-05-11 16:51:06.748997+00, null).

```

### setUpClass (tests.test_rls_policies_verification.RLSServiceRoleAccessTests)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 624, in setUpClass
    cls.test_data = cls.helper.create_test_data()
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 215, in create_test_data
    raise e
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 152, in create_test_data
    cursor.execute(f"""
    ~~~~~~~~~~~~~~^^^^^
        INSERT INTO team_members (team_id, user_id, role)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        VALUES ('{team_id}', '{user2_id}', 'member')
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    """)
    ^^^^
psycopg2.errors.CheckViolation: new row for relation "team_members" violates check constraint "team_members_role_check"
DETAIL:  Failing row contains (25eb2eb3-726f-4c00-8cfa-4bca839d524a, 33333333-3333-3333-3333-333333333333, 00000000-0000-0000-0000-000000000002, member, 2025-05-11 16:51:08.907413+00, 2025-05-11 16:51:08.907413+00, null).

```

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

```
