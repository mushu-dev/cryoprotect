'use client'

import { useState } from 'react'
import { useSession } from 'next-auth/react'
import { useRouter } from 'next/navigation'
import { Switch } from '@/components/ui/switch'
import { 
  Card, 
  CardContent, 
  CardDescription, 
  CardHeader, 
  CardTitle 
} from '@/components/ui/card'
import { 
  Tabs, 
  TabsContent, 
  TabsList, 
  TabsTrigger 
} from '@/components/ui/tabs'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import { 
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert'
import { Moon, Sun, Monitor, Bell, Globe, Trash } from 'lucide-react'

export default function SettingsPage() {
  const { data: session } = useSession()
  const router = useRouter()
  
  const [theme, setTheme] = useState('system')
  const [notifications, setNotifications] = useState({
    email: true,
    browser: true,
    updates: true,
    newsletter: false,
  })
  const [language, setLanguage] = useState('en')
  
  const [success, setSuccess] = useState<string | null>(null)
  const [error, setError] = useState<string | null>(null)
  
  // Redirect to login if not authenticated
  if (!session) {
    router.push('/auth/signin?callbackUrl=/settings')
    return null
  }
  
  // Handler for theme change
  const handleThemeChange = (value: string) => {
    setTheme(value)
    // In a real implementation, save this to user preferences in the backend
  }
  
  // Handler for notification toggle
  const handleNotificationToggle = (key: keyof typeof notifications) => {
    setNotifications({
      ...notifications,
      [key]: !notifications[key]
    })
  }
  
  // Handler for language change
  const handleLanguageChange = (value: string) => {
    setLanguage(value)
    // In a real implementation, save this to user preferences in the backend
  }
  
  // Handler for saving settings
  const handleSaveSettings = () => {
    try {
      // In a real implementation, save all settings to the backend
      setSuccess('Settings saved successfully')
      setError(null)
    } catch (err) {
      setError('Failed to save settings')
      setSuccess(null)
    }
  }
  
  return (
    <div className="max-w-4xl mx-auto">
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">Settings</h1>
        <p className="text-muted-foreground mt-1">
          Manage your application preferences and settings
        </p>
      </div>
      
      <Tabs defaultValue="appearance">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="appearance">Appearance</TabsTrigger>
          <TabsTrigger value="notifications">Notifications</TabsTrigger>
          <TabsTrigger value="account">Account</TabsTrigger>
        </TabsList>
        
        <TabsContent value="appearance" className="mt-6 space-y-4">
          {success && (
            <Alert>
              <AlertTitle>Success</AlertTitle>
              <AlertDescription>{success}</AlertDescription>
            </Alert>
          )}
          
          <Card>
            <CardHeader>
              <CardTitle>Appearance</CardTitle>
              <CardDescription>
                Customize the appearance of the application
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <div>
                  <h3 className="text-lg font-medium">Theme</h3>
                  <p className="text-sm text-muted-foreground">
                    Select a theme for the application
                  </p>
                </div>
                
                <div className="grid grid-cols-3 gap-4">
                  <Card 
                    className={`cursor-pointer ${theme === 'light' ? 'border-primary' : 'border-border'}`}
                    onClick={() => handleThemeChange('light')}
                  >
                    <CardContent className="flex flex-col items-center justify-center pt-6">
                      <Sun className="h-10 w-10 mb-2 text-yellow-500" />
                      <h4 className="font-medium">Light</h4>
                    </CardContent>
                  </Card>
                  
                  <Card
                    className={`cursor-pointer ${theme === 'dark' ? 'border-primary' : 'border-border'}`}
                    onClick={() => handleThemeChange('dark')}
                  >
                    <CardContent className="flex flex-col items-center justify-center pt-6">
                      <Moon className="h-10 w-10 mb-2 text-blue-500" />
                      <h4 className="font-medium">Dark</h4>
                    </CardContent>
                  </Card>
                  
                  <Card
                    className={`cursor-pointer ${theme === 'system' ? 'border-primary' : 'border-border'}`}
                    onClick={() => handleThemeChange('system')}
                  >
                    <CardContent className="flex flex-col items-center justify-center pt-6">
                      <Monitor className="h-10 w-10 mb-2 text-gray-500" />
                      <h4 className="font-medium">System</h4>
                    </CardContent>
                  </Card>
                </div>
              </div>
              
              <div className="space-y-4">
                <div>
                  <h3 className="text-lg font-medium">Language</h3>
                  <p className="text-sm text-muted-foreground">
                    Select your preferred language
                  </p>
                </div>
                
                <Select value={language} onValueChange={handleLanguageChange}>
                  <SelectTrigger className="w-[200px]">
                    <Globe className="h-4 w-4 mr-2" />
                    <SelectValue placeholder="Select language" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="en">English</SelectItem>
                    <SelectItem value="fr">Français</SelectItem>
                    <SelectItem value="de">Deutsch</SelectItem>
                    <SelectItem value="es">Español</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              
              <Button onClick={handleSaveSettings}>Save Changes</Button>
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="notifications" className="mt-6 space-y-4">
          {success && (
            <Alert>
              <AlertTitle>Success</AlertTitle>
              <AlertDescription>{success}</AlertDescription>
            </Alert>
          )}
          
          <Card>
            <CardHeader>
              <CardTitle>Notifications</CardTitle>
              <CardDescription>
                Configure how you receive notifications
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <h3 className="text-lg font-medium">Notification Settings</h3>
                
                <div className="space-y-4">
                  <div className="flex items-center justify-between">
                    <div className="space-y-0.5">
                      <Label htmlFor="email-notifications">Email Notifications</Label>
                      <p className="text-sm text-muted-foreground">
                        Receive notifications via email
                      </p>
                    </div>
                    <Switch
                      id="email-notifications"
                      checked={notifications.email}
                      onCheckedChange={() => handleNotificationToggle('email')}
                    />
                  </div>
                  
                  <div className="flex items-center justify-between">
                    <div className="space-y-0.5">
                      <Label htmlFor="browser-notifications">Browser Notifications</Label>
                      <p className="text-sm text-muted-foreground">
                        Receive notifications in the browser
                      </p>
                    </div>
                    <Switch
                      id="browser-notifications"
                      checked={notifications.browser}
                      onCheckedChange={() => handleNotificationToggle('browser')}
                    />
                  </div>
                  
                  <div className="flex items-center justify-between">
                    <div className="space-y-0.5">
                      <Label htmlFor="updates-notifications">Product Updates</Label>
                      <p className="text-sm text-muted-foreground">
                        Receive notifications about new features and updates
                      </p>
                    </div>
                    <Switch
                      id="updates-notifications"
                      checked={notifications.updates}
                      onCheckedChange={() => handleNotificationToggle('updates')}
                    />
                  </div>
                  
                  <div className="flex items-center justify-between">
                    <div className="space-y-0.5">
                      <Label htmlFor="newsletter-notifications">Newsletter</Label>
                      <p className="text-sm text-muted-foreground">
                        Receive our monthly newsletter
                      </p>
                    </div>
                    <Switch
                      id="newsletter-notifications"
                      checked={notifications.newsletter}
                      onCheckedChange={() => handleNotificationToggle('newsletter')}
                    />
                  </div>
                </div>
              </div>
              
              <Button onClick={handleSaveSettings}>Save Changes</Button>
            </CardContent>
          </Card>
        </TabsContent>
        
        <TabsContent value="account" className="mt-6 space-y-4">
          {error && (
            <Alert variant="destructive">
              <AlertTitle>Error</AlertTitle>
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}
          
          <Card>
            <CardHeader>
              <CardTitle>Account Settings</CardTitle>
              <CardDescription>
                Manage your account settings and preferences
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-4">
                <h3 className="text-lg font-medium">Account Information</h3>
                
                <div className="space-y-1">
                  <p className="text-sm font-medium">Email</p>
                  <p className="text-sm text-muted-foreground">{session.user?.email}</p>
                </div>
                
                <div className="space-y-1">
                  <p className="text-sm font-medium">Name</p>
                  <p className="text-sm text-muted-foreground">{session.user?.name || 'Not provided'}</p>
                </div>
                
                <div className="pt-2">
                  <Button variant="outline" asChild>
                    <Link href="/profile">Edit Profile</Link>
                  </Button>
                </div>
              </div>
              
              <div className="space-y-4 pt-4 border-t">
                <h3 className="text-lg font-medium text-destructive">Danger Zone</h3>
                
                <div className="space-y-2">
                  <h4 className="text-sm font-medium">Delete Account</h4>
                  <p className="text-sm text-muted-foreground">
                    Once you delete your account, there is no going back. Please be certain.
                  </p>
                  <Button variant="destructive">
                    <Trash className="h-4 w-4 mr-2" />
                    Delete Account
                  </Button>
                </div>
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  )
}

function Link({ href, children, ...props }: { href: string; children: React.ReactNode; [key: string]: any }) {
  return (
    <a href={href} {...props}>
      {children}
    </a>
  )
}