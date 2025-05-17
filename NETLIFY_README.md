# CryoProtect - Netlify Deployment

This project is configured to deploy on Netlify with the following settings:

## Configuration

The Netlify build should be configured as follows:

- **Base directory**: `frontend`
- **Build command**: `npm run build`
- **Publish directory**: `.next`
- **Node.js version**: 20.x or higher

## Environment Variables

The following environment variables should be set in Netlify:

- `NEXT_PUBLIC_API_URL`: `https://cryoprotect-8030e4025428.herokuapp.com/v1`
- `NEXT_PUBLIC_USE_MOCK_DATA`: `false`
- `NEXT_PUBLIC_ENABLE_API_LOGGING`: `true`
- `NEXT_PUBLIC_ENVIRONMENT`: `production`
- `NEXTAUTH_SECRET`: A random string for auth session encryption
- `NEXTAUTH_URL`: Your Netlify site URL
- `PROTECTION_BYPASS`: `TAt23KbtFE8dkZobJU3hpgTP4L5ja07V`
- `NEXT_PUBLIC_FRONTEND_PROTECTION_BYPASS`: `TAt23KbtFE8dkZobJU3hpgTP4L5ja07V`

## Important Notes

- The frontend code is in the `frontend` directory
- The site uses Next.js, not Svelte
- The legacy peer dependencies flag is enabled in the frontend's .npmrc file