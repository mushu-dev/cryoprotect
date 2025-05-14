release: python setup_database.py && python populate_sample_data.py
web: gunicorn simple_app:app --log-file -
