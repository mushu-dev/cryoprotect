import sys

modules_to_check = [
    'flask',
    'flask_restful',
    'flask_cors',
    'psutil',
    'pandas',
    'numpy',
    'matplotlib',
    'seaborn',
    'reportlab',
    'xlsxwriter',
    'requests',
    'sqlalchemy',
    'psycopg2',
    'jwt',
    'marshmallow',
    'rdkit',
    'rdkit.Chem'
]

for module_name in modules_to_check:
    try:
        print(f"Importing {module_name}...", end="")
        __import__(module_name)
        print(" ✓")
    except ImportError as e:
        print(f" ✗ Error: {e}")

print("\nAll imports completed.")
