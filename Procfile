release: python setup_database.py
web: python -c "import heroku_app_startup" && gunicorn app:app --log-file -
