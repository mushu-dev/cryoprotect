'use client';

import { useState } from 'react';
import { AnalyticsToggle } from '@/components/analytics/AnalyticsConsent';

export default function SettingsClientPage() {
  const [activeTab, setActiveTab] = useState('general');

  return (
    <div className="space-y-8">
      <div className="flex flex-col sm:flex-row gap-4">
        {/* Tabs */}
        <div className="w-full sm:w-64 flex flex-row sm:flex-col gap-2 p-2 border-b sm:border-b-0 sm:border-r border-gray-200 dark:border-gray-700">
          <button 
            className={`text-left px-3 py-2 rounded-md text-sm font-medium ${
              activeTab === 'general' 
                ? 'bg-gray-100 dark:bg-gray-800 text-gray-900 dark:text-white' 
                : 'text-gray-500 hover:text-gray-900 dark:hover:text-white hover:bg-gray-50 dark:hover:bg-gray-800'
            }`}
            onClick={() => setActiveTab('general')}
          >
            General
          </button>
          <button 
            className={`text-left px-3 py-2 rounded-md text-sm font-medium ${
              activeTab === 'appearance' 
                ? 'bg-gray-100 dark:bg-gray-800 text-gray-900 dark:text-white' 
                : 'text-gray-500 hover:text-gray-900 dark:hover:text-white hover:bg-gray-50 dark:hover:bg-gray-800'
            }`}
            onClick={() => setActiveTab('appearance')}
          >
            Appearance
          </button>
          <button 
            className={`text-left px-3 py-2 rounded-md text-sm font-medium ${
              activeTab === 'privacy' 
                ? 'bg-gray-100 dark:bg-gray-800 text-gray-900 dark:text-white' 
                : 'text-gray-500 hover:text-gray-900 dark:hover:text-white hover:bg-gray-50 dark:hover:bg-gray-800'
            }`}
            onClick={() => setActiveTab('privacy')}
          >
            Privacy
          </button>
          <button 
            className={`text-left px-3 py-2 rounded-md text-sm font-medium ${
              activeTab === 'notifications' 
                ? 'bg-gray-100 dark:bg-gray-800 text-gray-900 dark:text-white' 
                : 'text-gray-500 hover:text-gray-900 dark:hover:text-white hover:bg-gray-50 dark:hover:bg-gray-800'
            }`}
            onClick={() => setActiveTab('notifications')}
          >
            Notifications
          </button>
        </div>

        {/* Content */}
        <div className="flex-1 p-2">
          {activeTab === 'general' && (
            <div className="space-y-6">
              <h3 className="text-lg font-medium">General Settings</h3>
              <div className="border-b border-gray-200 dark:border-gray-700 pb-4">
                <h4 className="text-sm font-medium mb-1">Email Address</h4>
                <p className="text-sm text-gray-500 dark:text-gray-400 mb-2">
                  Update your email address
                </p>
                <div className="flex items-center gap-2">
                  <input 
                    type="email" 
                    className="flex h-9 rounded-md border border-input bg-background px-3 py-1 text-sm shadow-sm transition-colors file:border-0 file:bg-transparent file:text-sm file:font-medium placeholder:text-muted-foreground focus-visible:outline-none focus-visible:ring-1 focus-visible:ring-ring disabled:cursor-not-allowed disabled:opacity-50 flex-1"
                    placeholder="Enter your email" 
                  />
                  <button className="bg-blue-600 hover:bg-blue-700 text-white px-3 py-1 rounded-md text-sm">
                    Update
                  </button>
                </div>
              </div>
            </div>
          )}

          {activeTab === 'appearance' && (
            <div className="space-y-6">
              <h3 className="text-lg font-medium">Appearance Settings</h3>
              <div className="border-b border-gray-200 dark:border-gray-700 pb-4">
                <h4 className="text-sm font-medium mb-1">Theme</h4>
                <p className="text-sm text-gray-500 dark:text-gray-400 mb-2">
                  Choose your preferred theme
                </p>
                <div className="flex gap-2">
                  <button className="border border-gray-200 dark:border-gray-700 px-3 py-1 rounded-md text-sm">
                    Light
                  </button>
                  <button className="border border-gray-200 dark:border-gray-700 px-3 py-1 rounded-md text-sm">
                    Dark
                  </button>
                  <button className="border border-gray-200 dark:border-gray-700 px-3 py-1 rounded-md text-sm">
                    System
                  </button>
                </div>
              </div>
            </div>
          )}

          {activeTab === 'privacy' && (
            <div className="space-y-6">
              <h3 className="text-lg font-medium">Privacy Settings</h3>
              
              <div className="border-b border-gray-200 dark:border-gray-700 pb-4">
                <AnalyticsToggle />
              </div>
              
              <div className="border-b border-gray-200 dark:border-gray-700 pb-4">
                <h4 className="text-sm font-medium mb-1">Data Sharing</h4>
                <p className="text-sm text-gray-500 dark:text-gray-400 mb-2">
                  Control how your data is shared
                </p>
                <div className="flex items-center justify-between py-3">
                  <div>
                    <h3 className="text-sm font-medium">Share usage statistics</h3>
                    <p className="text-sm text-gray-500 dark:text-gray-400">
                      Allow anonymous sharing of how you use the application
                    </p>
                  </div>
                  <button
                    type="button"
                    role="switch"
                    aria-checked={false}
                    className="relative inline-flex h-6 w-11 items-center rounded-full transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 bg-gray-200 dark:bg-gray-700"
                  >
                    <span className="sr-only">Enable data sharing</span>
                    <span
                      className="inline-block h-4 w-4 transform rounded-full bg-white transition-transform translate-x-1"
                    />
                  </button>
                </div>
              </div>
            </div>
          )}

          {activeTab === 'notifications' && (
            <div className="space-y-6">
              <h3 className="text-lg font-medium">Notification Settings</h3>
              <div className="border-b border-gray-200 dark:border-gray-700 pb-4">
                <h4 className="text-sm font-medium mb-1">Email Notifications</h4>
                <p className="text-sm text-gray-500 dark:text-gray-400 mb-2">
                  Manage your email notification preferences
                </p>
                <div className="space-y-2">
                  <div className="flex items-center">
                    <input 
                      type="checkbox" 
                      id="updates" 
                      className="h-4 w-4 rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                    />
                    <label htmlFor="updates" className="ml-2 text-sm text-gray-700 dark:text-gray-300">
                      Product updates and announcements
                    </label>
                  </div>
                  <div className="flex items-center">
                    <input 
                      type="checkbox" 
                      id="security" 
                      className="h-4 w-4 rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                    />
                    <label htmlFor="security" className="ml-2 text-sm text-gray-700 dark:text-gray-300">
                      Security alerts
                    </label>
                  </div>
                </div>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}