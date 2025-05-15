# Connecting Heroku to RDKit Service

Now that your RDKit service is deployed at https://cryoprotect-rdkit.fly.dev/, follow these steps to connect it to your Heroku application:

## 1. Set Environment Variables in Heroku

```bash
# Set the RDKit service URL in your Heroku app
heroku config:set RDKIT_SERVICE_URL=https://cryoprotect-rdkit.fly.dev -a your-heroku-app-name

# Optional: Set an API key for security
API_KEY=$(openssl rand -hex 16)
echo "Generated API key: $API_KEY"

# Set the API key in both applications
flyctl secrets set API_KEY=$API_KEY
heroku config:set RDKIT_API_KEY=$API_KEY -a your-heroku-app-name
```

## 2. Add the RDKit Client to Your Heroku App

Copy the `heroku-client.py` file to your Heroku application:

```bash
# Copy from the current location to your application directory
cp /home/mushu/Projects/CryoProtect/heroku-client.py /path/to/your/heroku/app/rdkit_client.py
```

## 3. Use the Client in Your Heroku App

Update your application code to use the RDKit client:

```python
from rdkit_client import RDKitServiceClient

# Initialize the client
rdkit_client = RDKitServiceClient()

# Example: Calculate properties for a molecule
@app.route('/api/molecules/<int:molecule_id>/properties')
def get_molecule_properties(molecule_id):
    # Get the molecule from the database
    molecule = get_molecule_by_id(molecule_id)
    
    # Calculate properties using the RDKit service
    properties = rdkit_client.calculate_properties(molecule.smiles)
    
    return jsonify(properties)
```

## 4. Test the Integration

1. Deploy your updated Heroku app:
   ```bash
   git add .
   git commit -m "Integrate RDKit service"
   git push heroku master
   ```

2. Test an API endpoint that uses the RDKit service:
   ```bash
   curl https://your-heroku-app.herokuapp.com/api/molecules/1/properties
   ```

## 5. Monitoring and Troubleshooting

### Monitor the RDKit service:
```bash
flyctl logs
```

### Check if the RDKit service is running:
```bash
curl https://cryoprotect-rdkit.fly.dev/health
```

### Test RDKit functionality directly:
```bash
curl https://cryoprotect-rdkit.fly.dev/rdkit-check
```

### Test a specific endpoint:
```bash
curl -X POST https://cryoprotect-rdkit.fly.dev/api/calculate-properties \
  -H "Content-Type: application/json" \
  -d '{"molecule_data": "CCO"}'
```

## 6. Scaling Considerations

- If your application needs more RDKit processing power, you can scale the Fly.io service:
  ```bash
  flyctl scale count 2
  ```

- For cost optimization, the service is configured to auto-stop when not in use and auto-start when needed.

## 7. Backup Plan

If the RDKit service becomes unavailable, the client will automatically fall back to mock implementations, ensuring your application remains functional with limited capabilities.

## 8. Security Notes

- The service is secured with HTTPS
- Add API key authentication for production use
- Consider IP restrictions for additional security