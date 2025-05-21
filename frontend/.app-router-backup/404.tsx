export default function NotFound() {
  return (
    <div className="flex flex-col items-center justify-center min-h-[70vh]">
      <h1 className="text-4xl font-bold">404 - Page Not Found</h1>
      <p className="mt-4 text-lg text-muted-foreground">
        The page you are looking for does not exist.
      </p>
      <a 
        href="/"
        className="mt-6 px-4 py-2 rounded bg-primary text-primary-foreground hover:bg-primary/90"
      >
        Return Home
      </a>
    </div>
  )
}
