'use client';

import { useState, useEffect } from 'react'
import { useSession } from 'next-auth/react'
import { type Session } from 'next-auth'
import { useRouter } from 'next/navigation'
import { z } from 'zod'
import { useForm } from 'react-hook-form'
import { zodResolver } from '@hookform/resolvers/zod'
import { 
  Card, 
  CardContent, 
  CardDescription, 
  CardFooter, 
  CardHeader, 
  CardTitle 
} from '@/components/ui/card'
import { 
  Tabs, 
  TabsContent, 
  TabsList, 
  TabsTrigger 
} from '@/components/ui/tabs'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import { Textarea } from '@/components/ui/textarea'
import { Avatar, AvatarFallback, AvatarImage } from '@/components/ui/avatar'
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert'
import { User, Mail, Building, Edit, Save, Shield, Key, CheckCircle2 } from 'lucide-react'

// Custom session interface to include our additional properties
interface CustomSession extends Session {
  accessToken?: string;
  error?: string;
}

// Validation schema for profile update
const profileSchema = z.object({
  name: z.string().min(2, 'Name must be at least 2 characters'),
  organization: z.string().optional(),
  bio: z.string().optional(),
})

type ProfileFormValues = z.infer<typeof profileSchema>

// Validation schema for password change
const passwordSchema = z.object({
  currentPassword: z.string().min(1, 'Current password is required'),
  newPassword: z.string().min(8, 'New password must be at least 8 characters'),
  confirmPassword: z.string().min(8, 'Password confirmation is required'),
}).refine((data) => data.newPassword === data.confirmPassword, {
  message: "Passwords don't match",
  path: ["confirmPassword"],
})

type PasswordFormValues = z.infer<typeof passwordSchema>

export default function ProfileClient() {
  const { data: session, update } = useSession()
  const router = useRouter()
  
  const [isUpdating, setIsUpdating] = useState(false)
  const [isChangingPassword, setIsChangingPassword] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [success, setSuccess] = useState<string | null>(null)
  
  // Redirect to login if not authenticated
  if (!session) {
    router.push('/auth/signin?callbackUrl=/profile')
    return null
  }
  
  // Cast session to our custom type
  const customSession = session as CustomSession
  
  // Profile form setup
  const profileForm = useForm<ProfileFormValues>({
    resolver: zodResolver(profileSchema),
    defaultValues: {
      name: session.user?.name || '',
      organization: '',
      bio: '',
    },
  })
  
  // Password form setup
  const passwordForm = useForm<PasswordFormValues>({
    resolver: zodResolver(passwordSchema),
    defaultValues: {
      currentPassword: '',
      newPassword: '',
      confirmPassword: '',
    },
  })
  
  // Handle profile update
  const onUpdateProfile = async (values: ProfileFormValues) => {
    try {
      setIsUpdating(true)
      setError(null)
      setSuccess(null)
      
      // In a real implementation, make an API request to update the user's profile
      const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL}/users/profile`, {
        method: 'PATCH',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${customSession.accessToken}`,
        },
        body: JSON.stringify(values),
      })
      
      const data = await response.json()
      
      if (!response.ok) {
        throw new Error(data.message || 'Failed to update profile')
      }
      
      // Update the session with the new user data
      await update({
        ...session,
        user: {
          ...session.user,
          name: values.name,
        },
      })
      
      setSuccess('Your profile has been updated successfully!')
    } catch (error) {
      console.error('Error updating profile:', error)
      setError(error instanceof Error ? error.message : 'An error occurred while updating your profile')
    } finally {
      setIsUpdating(false)
    }
  }
  
  // Handle password change
  const onChangePassword = async (values: PasswordFormValues) => {
    try {
      setIsChangingPassword(true)
      setError(null)
      setSuccess(null)
      
      // In a real implementation, make an API request to change the user's password
      const response = await fetch(`${process.env.NEXT_PUBLIC_API_URL}/users/password`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${customSession.accessToken}`,
        },
        body: JSON.stringify({
          currentPassword: values.currentPassword,
          newPassword: values.newPassword,
        }),
      })
      
      const data = await response.json()
      
      if (!response.ok) {
        throw new Error(data.message || 'Failed to change password')
      }
      
      // Reset the form
      passwordForm.reset()
      
      setSuccess('Your password has been changed successfully. Please use your new password next time you log in.')
    } catch (error) {
      console.error('Error changing password:', error)
      setError(error instanceof Error ? error.message : 'An error occurred while changing your password')
    } finally {
      setIsChangingPassword(false)
    }
  }
  
  // Generate avatar fallback from user's name
  const getInitials = (name: string) => {
    const names = name.split(' ')
    if (names.length === 1) return names[0].charAt(0).toUpperCase()
    return (names[0].charAt(0) + names[names.length - 1].charAt(0)).toUpperCase()
  }
  
  const userInitials = session.user?.name ? getInitials(session.user.name) : 'U'
  
  // Effect for handling success messages - auto-hide after 5 seconds
  useEffect(() => {
    if (success) {
      const timer = setTimeout(() => {
        setSuccess(null);
      }, 5000);
      return () => clearTimeout(timer);
    }
  }, [success]);

  return (
    <div className="container mx-auto px-4 py-6 max-w-5xl">
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">Profile</h1>
        <p className="text-muted-foreground mt-1">
          Manage your account settings and preferences
        </p>
      </div>
      
      <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
        <div className="lg:col-span-1">
          <Card className="sticky top-20">
            <CardContent className="pt-6">
              <div className="flex flex-col items-center space-y-4">
                <Avatar className="h-24 w-24 border-2 border-primary/10">
                  <AvatarImage src={session.user?.image || ''} alt={session.user?.name || 'User'} />
                  <AvatarFallback className="text-xl bg-primary/5">{userInitials}</AvatarFallback>
                </Avatar>
                
                <div className="text-center">
                  <h2 className="text-xl font-semibold">{session.user?.name}</h2>
                  <p className="text-sm text-muted-foreground mt-1">
                    {session.user?.email}
                  </p>
                </div>
                
                <Button variant="outline" className="w-full">
                  <Edit className="h-4 w-4 mr-2" />
                  Change Photo
                </Button>
              </div>
            </CardContent>
          </Card>
        </div>
        
        <div className="lg:col-span-3">
          <Tabs defaultValue="profile" className="w-full">
            <TabsList className="w-full sm:w-auto flex sm:inline-flex overflow-x-auto">
              <TabsTrigger value="profile" className="flex-1 sm:flex-initial">
                <User className="h-4 w-4 mr-2" />
                Profile
              </TabsTrigger>
              <TabsTrigger value="security" className="flex-1 sm:flex-initial">
                <Shield className="h-4 w-4 mr-2" />
                Security
              </TabsTrigger>
            </TabsList>
            
            <TabsContent value="profile" className="mt-6 space-y-6">
              {error && (
                <Alert variant="destructive" className="animate-in fade-in-50">
                  <AlertTitle>Error</AlertTitle>
                  <AlertDescription>{error}</AlertDescription>
                </Alert>
              )}
              
              {success && (
                <Alert className="bg-green-50 dark:bg-green-900/20 border-green-200 dark:border-green-900/30 animate-in fade-in-50">
                  <CheckCircle2 className="h-4 w-4 text-green-600 dark:text-green-400" />
                  <AlertTitle className="text-green-800 dark:text-green-200">Success</AlertTitle>
                  <AlertDescription className="text-green-700 dark:text-green-300">{success}</AlertDescription>
                </Alert>
              )}
              
              <Card className="shadow-sm border-primary/10">
                <CardHeader className="pb-4">
                  <CardTitle>Personal Information</CardTitle>
                  <CardDescription>
                    Update your personal information and profile details
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <form onSubmit={profileForm.handleSubmit(onUpdateProfile)} className="space-y-6">
                    <div className="space-y-3">
                      <div className="flex items-center">
                        <User className="h-4 w-4 mr-2 text-muted-foreground" />
                        <Label htmlFor="name">Name</Label>
                      </div>
                      <Input
                        id="name"
                        placeholder="Your full name"
                        className="transition-all focus-visible:ring-primary/50"
                        {...profileForm.register('name')}
                      />
                      {profileForm.formState.errors.name && (
                        <p className="text-destructive text-sm mt-1 animate-in slide-in-from-left-1">
                          {profileForm.formState.errors.name.message}
                        </p>
                      )}
                    </div>
                    
                    <div className="space-y-3">
                      <div className="flex items-center">
                        <Mail className="h-4 w-4 mr-2 text-muted-foreground" />
                        <Label htmlFor="email">Email</Label>
                      </div>
                      <Input
                        id="email"
                        value={session.user?.email || ''}
                        disabled
                        readOnly
                        className="opacity-70"
                      />
                      <p className="text-xs text-muted-foreground mt-1">
                        Email address cannot be changed
                      </p>
                    </div>
                    
                    <div className="space-y-3">
                      <div className="flex items-center">
                        <Building className="h-4 w-4 mr-2 text-muted-foreground" />
                        <Label htmlFor="organization">Organization</Label>
                      </div>
                      <Input
                        id="organization"
                        placeholder="Your organization"
                        className="transition-all focus-visible:ring-primary/50"
                        {...profileForm.register('organization')}
                      />
                    </div>
                    
                    <div className="space-y-3">
                      <Label htmlFor="bio" className="flex items-start gap-2">
                        <span>Bio</span>
                        <span className="text-xs text-muted-foreground">(optional)</span>
                      </Label>
                      <Textarea
                        id="bio"
                        placeholder="Tell us about yourself"
                        rows={4}
                        className="transition-all focus-visible:ring-primary/50 min-h-[100px] resize-y"
                        {...profileForm.register('bio')}
                      />
                    </div>
                    
                    <div className="pt-2">
                      <Button 
                        type="submit" 
                        disabled={isUpdating}
                        className="w-full sm:w-auto transition-all"
                      >
                        {isUpdating ? (
                          <>
                            <svg className="animate-spin -ml-1 mr-2 h-4 w-4 text-white" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                              <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                              <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                            </svg>
                            Saving...
                          </>
                        ) : (
                          <>
                            <Save className="h-4 w-4 mr-2" />
                            Save Changes
                          </>
                        )}
                      </Button>
                    </div>
                  </form>
                </CardContent>
              </Card>
            </TabsContent>
            
            <TabsContent value="security" className="mt-6 space-y-6">
              {error && (
                <Alert variant="destructive" className="animate-in fade-in-50">
                  <AlertTitle>Error</AlertTitle>
                  <AlertDescription>{error}</AlertDescription>
                </Alert>
              )}
              
              {success && (
                <Alert className="bg-green-50 dark:bg-green-900/20 border-green-200 dark:border-green-900/30 animate-in fade-in-50">
                  <CheckCircle2 className="h-4 w-4 text-green-600 dark:text-green-400" />
                  <AlertTitle className="text-green-800 dark:text-green-200">Success</AlertTitle>
                  <AlertDescription className="text-green-700 dark:text-green-300">{success}</AlertDescription>
                </Alert>
              )}
              
              <Card className="shadow-sm border-primary/10">
                <CardHeader className="pb-4">
                  <CardTitle>Password</CardTitle>
                  <CardDescription>
                    Change your password to keep your account secure
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <form onSubmit={passwordForm.handleSubmit(onChangePassword)} className="space-y-6">
                    <div className="space-y-3">
                      <div className="flex items-center">
                        <Key className="h-4 w-4 mr-2 text-muted-foreground" />
                        <Label htmlFor="currentPassword">Current Password</Label>
                      </div>
                      <Input
                        id="currentPassword"
                        type="password"
                        className="transition-all focus-visible:ring-primary/50"
                        {...passwordForm.register('currentPassword')}
                      />
                      {passwordForm.formState.errors.currentPassword && (
                        <p className="text-destructive text-sm mt-1 animate-in slide-in-from-left-1">
                          {passwordForm.formState.errors.currentPassword.message}
                        </p>
                      )}
                    </div>
                    
                    <div className="space-y-3">
                      <Label htmlFor="newPassword">New Password</Label>
                      <Input
                        id="newPassword"
                        type="password"
                        className="transition-all focus-visible:ring-primary/50"
                        {...passwordForm.register('newPassword')}
                      />
                      {passwordForm.formState.errors.newPassword && (
                        <p className="text-destructive text-sm mt-1 animate-in slide-in-from-left-1">
                          {passwordForm.formState.errors.newPassword.message}
                        </p>
                      )}
                    </div>
                    
                    <div className="space-y-3">
                      <Label htmlFor="confirmPassword">Confirm New Password</Label>
                      <Input
                        id="confirmPassword"
                        type="password"
                        className="transition-all focus-visible:ring-primary/50"
                        {...passwordForm.register('confirmPassword')}
                      />
                      {passwordForm.formState.errors.confirmPassword && (
                        <p className="text-destructive text-sm mt-1 animate-in slide-in-from-left-1">
                          {passwordForm.formState.errors.confirmPassword.message}
                        </p>
                      )}
                    </div>
                    
                    <div className="pt-2">
                      <Button 
                        type="submit" 
                        disabled={isChangingPassword}
                        className="w-full sm:w-auto transition-all"
                      >
                        {isChangingPassword ? (
                          <>
                            <svg className="animate-spin -ml-1 mr-2 h-4 w-4 text-white" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                              <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                              <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                            </svg>
                            Changing...
                          </>
                        ) : 'Change Password'}
                      </Button>
                    </div>
                  </form>
                </CardContent>
              </Card>
              
              <Card className="shadow-sm border-primary/10 mt-8">
                <CardHeader className="pb-4">
                  <CardTitle>Two-Factor Authentication</CardTitle>
                  <CardDescription>
                    Add an extra layer of security to your account
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <p className="text-sm text-muted-foreground mb-6">
                    Two-factor authentication adds an extra layer of security to your account by requiring more than just a password to sign in.
                    When enabled, you'll need to provide a verification code from your phone when logging in.
                  </p>
                  <Button variant="outline" className="transition-all hover:bg-primary/5">
                    Enable Two-Factor Authentication
                  </Button>
                </CardContent>
              </Card>
            </TabsContent>
          </Tabs>
        </div>
      </div>
    </div>
  )
}