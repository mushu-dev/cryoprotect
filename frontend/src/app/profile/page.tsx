'use client'

import { useState } from 'react'
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
import { User, Mail, Building, Edit, Save, Shield, Key } from 'lucide-react'

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

export default function ProfilePage() {
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
      
      setSuccess('Profile updated successfully')
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
      
      setSuccess('Password changed successfully')
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
  
  return (
    <div className="max-w-4xl mx-auto">
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">Profile</h1>
        <p className="text-muted-foreground mt-1">
          Manage your account settings and preferences
        </p>
      </div>
      
      <div className="grid grid-cols-1 md:grid-cols-4 gap-6">
        <div className="md:col-span-1">
          <Card>
            <CardContent className="pt-6">
              <div className="flex flex-col items-center space-y-4">
                <Avatar className="h-24 w-24">
                  <AvatarImage src={session.user?.image || ''} alt={session.user?.name || 'User'} />
                  <AvatarFallback className="text-xl">{userInitials}</AvatarFallback>
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
        
        <div className="md:col-span-3">
          <Tabs defaultValue="profile">
            <TabsList className="grid w-full grid-cols-2">
              <TabsTrigger value="profile">
                <User className="h-4 w-4 mr-2" />
                Profile
              </TabsTrigger>
              <TabsTrigger value="security">
                <Shield className="h-4 w-4 mr-2" />
                Security
              </TabsTrigger>
            </TabsList>
            
            <TabsContent value="profile" className="mt-6 space-y-4">
              {error && (
                <Alert variant="destructive">
                  <AlertTitle>Error</AlertTitle>
                  <AlertDescription>{error}</AlertDescription>
                </Alert>
              )}
              
              {success && (
                <Alert>
                  <AlertTitle>Success</AlertTitle>
                  <AlertDescription>{success}</AlertDescription>
                </Alert>
              )}
              
              <Card>
                <CardHeader>
                  <CardTitle>Personal Information</CardTitle>
                  <CardDescription>
                    Update your personal information and profile details
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <form onSubmit={profileForm.handleSubmit(onUpdateProfile)} className="space-y-4">
                    <div className="space-y-2">
                      <div className="flex items-center">
                        <User className="h-4 w-4 mr-2 text-muted-foreground" />
                        <Label htmlFor="name">Name</Label>
                      </div>
                      <Input
                        id="name"
                        placeholder="Your full name"
                        {...profileForm.register('name')}
                      />
                      {profileForm.formState.errors.name && (
                        <p className="text-destructive text-sm">
                          {profileForm.formState.errors.name.message}
                        </p>
                      )}
                    </div>
                    
                    <div className="space-y-2">
                      <div className="flex items-center">
                        <Mail className="h-4 w-4 mr-2 text-muted-foreground" />
                        <Label htmlFor="email">Email</Label>
                      </div>
                      <Input
                        id="email"
                        value={session.user?.email || ''}
                        disabled
                        readOnly
                      />
                      <p className="text-xs text-muted-foreground">
                        Email address cannot be changed
                      </p>
                    </div>
                    
                    <div className="space-y-2">
                      <div className="flex items-center">
                        <Building className="h-4 w-4 mr-2 text-muted-foreground" />
                        <Label htmlFor="organization">Organization</Label>
                      </div>
                      <Input
                        id="organization"
                        placeholder="Your organization"
                        {...profileForm.register('organization')}
                      />
                    </div>
                    
                    <div className="space-y-2">
                      <Label htmlFor="bio">Bio</Label>
                      <Textarea
                        id="bio"
                        placeholder="Tell us about yourself"
                        rows={4}
                        {...profileForm.register('bio')}
                      />
                    </div>
                    
                    <Button type="submit" disabled={isUpdating}>
                      {isUpdating ? 'Saving...' : (
                        <>
                          <Save className="h-4 w-4 mr-2" />
                          Save Changes
                        </>
                      )}
                    </Button>
                  </form>
                </CardContent>
              </Card>
            </TabsContent>
            
            <TabsContent value="security" className="mt-6 space-y-4">
              {error && (
                <Alert variant="destructive">
                  <AlertTitle>Error</AlertTitle>
                  <AlertDescription>{error}</AlertDescription>
                </Alert>
              )}
              
              {success && (
                <Alert>
                  <AlertTitle>Success</AlertTitle>
                  <AlertDescription>{success}</AlertDescription>
                </Alert>
              )}
              
              <Card>
                <CardHeader>
                  <CardTitle>Password</CardTitle>
                  <CardDescription>
                    Change your password to keep your account secure
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <form onSubmit={passwordForm.handleSubmit(onChangePassword)} className="space-y-4">
                    <div className="space-y-2">
                      <div className="flex items-center">
                        <Key className="h-4 w-4 mr-2 text-muted-foreground" />
                        <Label htmlFor="currentPassword">Current Password</Label>
                      </div>
                      <Input
                        id="currentPassword"
                        type="password"
                        {...passwordForm.register('currentPassword')}
                      />
                      {passwordForm.formState.errors.currentPassword && (
                        <p className="text-destructive text-sm">
                          {passwordForm.formState.errors.currentPassword.message}
                        </p>
                      )}
                    </div>
                    
                    <div className="space-y-2">
                      <Label htmlFor="newPassword">New Password</Label>
                      <Input
                        id="newPassword"
                        type="password"
                        {...passwordForm.register('newPassword')}
                      />
                      {passwordForm.formState.errors.newPassword && (
                        <p className="text-destructive text-sm">
                          {passwordForm.formState.errors.newPassword.message}
                        </p>
                      )}
                    </div>
                    
                    <div className="space-y-2">
                      <Label htmlFor="confirmPassword">Confirm New Password</Label>
                      <Input
                        id="confirmPassword"
                        type="password"
                        {...passwordForm.register('confirmPassword')}
                      />
                      {passwordForm.formState.errors.confirmPassword && (
                        <p className="text-destructive text-sm">
                          {passwordForm.formState.errors.confirmPassword.message}
                        </p>
                      )}
                    </div>
                    
                    <Button type="submit" disabled={isChangingPassword}>
                      {isChangingPassword ? 'Changing...' : 'Change Password'}
                    </Button>
                  </form>
                </CardContent>
              </Card>
              
              <Card>
                <CardHeader>
                  <CardTitle>Two-Factor Authentication</CardTitle>
                  <CardDescription>
                    Add an extra layer of security to your account
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <p className="text-sm text-muted-foreground mb-4">
                    Two-factor authentication adds an extra layer of security to your account by requiring more than just a password to sign in.
                  </p>
                  <Button variant="outline">
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