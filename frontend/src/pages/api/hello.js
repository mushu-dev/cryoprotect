export default function handler(req, res) {
  res.status(200).json({
    name: 'CryoProtect Minimal API',
    message: 'API is working!',
    timestamp: new Date().toISOString(),
    environment: process.env.NEXT_PUBLIC_ENVIRONMENT || 'unknown'
  });
}