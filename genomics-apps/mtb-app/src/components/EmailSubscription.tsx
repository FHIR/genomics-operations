'use client';

import React, { useState } from 'react';
import { Mail, Check, X } from 'lucide-react';

export default function EmailSubscription() {
  const [email, setEmail] = useState('');
  const [status, setStatus] = useState<'idle' | 'submitted' | 'error'>('idle');
  const [message, setMessage] = useState('');

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    
    if (!email.trim()) {
      setStatus('error');
      setMessage('Please enter an email address');
      return;
    }

    // Basic email validation
    const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
    if (!emailRegex.test(email)) {
      setStatus('error');
      setMessage('Please enter a valid email address');
      return;
    }

    // Simulate submission (since no backend)
    console.log('Email subscription submitted:', email);
    setStatus('submitted');
    setMessage("Thank you! We'll keep you updated.");
    setEmail('');

    setTimeout(() => {
      setStatus('idle');
      setMessage('');
    }, 3000);
  };
  return (
    <div className="bg-gray-50 p-6 rounded-lg border border-gray-200">
      <div className="flex items-center gap-2 mb-4">
        <Mail className="text-blue-600" size={20} />
        <h3 className="text-lg font-semibold text-gray-800">Stay Updated</h3>
      </div>
      
      <p className="text-gray-600 mb-4">
        Get notified about new features and updates to the Molecular Tumor Board platform.
      </p>

      <form onSubmit={handleSubmit} className="space-y-3">
        <div>
          <input
            type="email"
            value={email}
            onChange={(e) => setEmail(e.target.value)}
            placeholder="Enter your email address"
            className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
            disabled={status === 'submitted'}
          />
        </div>

        <button
          type="submit"
          disabled={status === 'submitted'}
          className="w-full bg-blue-600 text-white py-2 px-4 rounded-md hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
        >
          {status === 'submitted' ? 'Subscribed!' : 'Subscribe for Updates'}
        </button>

        {message && (
          <div className={`flex items-center gap-2 text-sm ${
            status === 'error' ? 'text-red-600' : 'text-green-600'
          }`}>
            {status === 'error' ? <X size={16} /> : <Check size={16} />}
            {message}
          </div>
        )}
      </form>
    </div>
  );
}