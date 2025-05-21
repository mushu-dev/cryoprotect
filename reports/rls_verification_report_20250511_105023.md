# RLS Verification Report

Generated on: 2025-05-11 10:50:23

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
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 44, in create_test_data
    cursor.execute(f"""
    ~~~~~~~~~~~~~~^^^^^
        INSERT INTO auth.users (id, email, created_at, updated_at, role, aud, confirmed_at)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ...<8 lines>...
        )
        ^
    """)
    ^^^^
psycopg2.errors.GeneratedAlways: cannot insert a non-DEFAULT value into column "confirmed_at"
DETAIL:  Column "confirmed_at" is a generated column.

```

### setUpClass (tests.test_rls_policies_verification.RLSTableAccessPolicyTests)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 231, in setUpClass
    cls.test_data = cls.helper.create_test_data()
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 215, in create_test_data
    raise e
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 44, in create_test_data
    cursor.execute(f"""
    ~~~~~~~~~~~~~~^^^^^
        INSERT INTO auth.users (id, email, created_at, updated_at, role, aud, confirmed_at)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ...<8 lines>...
        )
        ^
    """)
    ^^^^
psycopg2.errors.GeneratedAlways: cannot insert a non-DEFAULT value into column "confirmed_at"
DETAIL:  Column "confirmed_at" is a generated column.

```

### setUpClass (tests.test_rls_policies_verification.RLSRelationshipPolicyTests)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 380, in setUpClass
    cls.test_data = cls.helper.create_test_data()
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 215, in create_test_data
    raise e
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 44, in create_test_data
    cursor.execute(f"""
    ~~~~~~~~~~~~~~^^^^^
        INSERT INTO auth.users (id, email, created_at, updated_at, role, aud, confirmed_at)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ...<8 lines>...
        )
        ^
    """)
    ^^^^
psycopg2.errors.GeneratedAlways: cannot insert a non-DEFAULT value into column "confirmed_at"
DETAIL:  Column "confirmed_at" is a generated column.

```

### setUpClass (tests.test_rls_policies_verification.RLSServiceRoleAccessTests)

```
Traceback (most recent call last):
  File "/home/mushu/Projects/CryoProtect/tests/test_rls_policies_verification.py", line 624, in setUpClass
    cls.test_data = cls.helper.create_test_data()
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 215, in create_test_data
    raise e
  File "/home/mushu/Projects/CryoProtect/tests/rls_test_helper.py", line 44, in create_test_data
    cursor.execute(f"""
    ~~~~~~~~~~~~~~^^^^^
        INSERT INTO auth.users (id, email, created_at, updated_at, role, aud, confirmed_at)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ...<8 lines>...
        )
        ^
    """)
    ^^^^
psycopg2.errors.GeneratedAlways: cannot insert a non-DEFAULT value into column "confirmed_at"
DETAIL:  Column "confirmed_at" is a generated column.

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

```
