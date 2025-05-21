import importlib
import sys

def check_dependency(name):
    try:
        importlib.import_module(name)
        return True
    except ImportError:
        return False

dependencies = [
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
    'pyjwt',
    'marshmallow',
    'rdkit'
]

print("Dependency Check Results:")
print("-------------------------")

missing = []
for dep in dependencies:
    if check_dependency(dep):
        print(f"✅ {dep} - Installed")
    else:
        print(f"❌ {dep} - Missing")
        missing.append(dep)

if missing:
    print("\nMissing Dependencies:")
    print("--------------------")
    for dep in missing:
        print(f"pip install {dep}")
else:
    print("\nAll dependencies are installed\!")
