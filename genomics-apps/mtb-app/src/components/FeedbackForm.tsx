'use client';

import React, { useState } from 'react';
import { MessageSquare, Send, Check, X } from 'lucide-react';

export default function FeedbackForm() {
  const [feedback, setFeedback] = useState('');
  const [email, setEmail] = useState('');
  const [category, setCategory] = useState('general');
  const [status, setStatus] = useState<'idle' | 'submitted' | 'error'>('idle');
  const [message, setMessage] = useState('');

  const categories = [
    { value: 'general', label: 'General Feedback' },
    { value: 'bug', label: 'Bug Report' },
    { value: 'feature', label: 'Feature Request' },
    { value: 'usability', label: 'Usability Issue' },
    { value: 'data', label: 'Data Quality' }
  ];

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    
    if (!feedback.trim()) {
      setStatus('error');
      setMessage('Please enter your feedback');
      return;
    }

    if (feedback.trim().length < 10) {
      setStatus('error');
      setMessage('Please provide more detailed feedback (at least 10 characters)');
      return;
    }

    // Basic email validation if provided
    if (email.trim() && !/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(email)) {
      setStatus('error');
      setMessage('Please enter a valid email address');
      return;
    }

    // Simulate submission (since no backend)
    const feedbackData = {
      feedback: feedback.trim(),
      email: email.trim(),
      category,
      timestamp: new Date().toISOString()
    };
    
    console.log('Feedback submitted:', feedbackData);
    setStatus('submitted');
    setMessage('Thank you for your feedback! We appreciate your input.');
    
    // Reset form
    setFeedback('');
    setEmail('');
    setCategory('general');

    setTimeout(() => {
      setStatus('idle');
      setMessage('');
    }, 4000);
  };

  return (
    <div className="bg-gray-50 p-6 rounded-lg border border-gray-200">
      <div className="flex items-center gap-2 mb-4">
        <MessageSquare className="text-green-600" size={20} />
        <h3 className="text-lg font-semibold text-gray-800">Share Your Feedback</h3>
      </div>
      
      <p className="text-gray-600 mb-4">
        Help us improve the platform by sharing your thoughts, reporting bugs, or suggesting new features.
      </p>

      <form onSubmit={handleSubmit} className="space-y-4">
        <div>
          <label htmlFor="category" className="block text-sm font-medium text-gray-700 mb-1">
            Category
          </label>
          <select
            id="category"
            value={category}
            onChange={(e) => setCategory(e.target.value)}
            className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-green-500 focus:border-transparent"
            disabled={status === 'submitted'}
          >
            {categories.map((cat) => (
              <option key={cat.value} value={cat.value}>
                {cat.label}
              </option>
            ))}
          </select>
        </div>

        <div>
          <label htmlFor="feedback" className="block text-sm font-medium text-gray-700 mb-1">
            Your Feedback <span className="text-red-500">*</span>
          </label>
          <textarea
            id="feedback"
            value={feedback}
            onChange={(e) => setFeedback(e.target.value)}
            placeholder="Tell us about your experience, report a bug, or suggest an improvement..."
            rows={4}
            className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-green-500 focus:border-transparent resize-vertical"
            disabled={status === 'submitted'}
          />
          <p className="text-xs text-gray-500 mt-1">
            {feedback.length}/500 characters
          </p>
        </div>

        <div>
          <label htmlFor="email" className="block text-sm font-medium text-gray-700 mb-1">
            Email Address (optional)
          </label>
          <input
            type="email"
            id="email"
            value={email}
            onChange={(e) => setEmail(e.target.value)}
            placeholder="your.email@example.com"
            className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-green-500 focus:border-transparent"
            disabled={status === 'submitted'}
          />
          <p className="text-xs text-gray-500 mt-1">
            Provide your email if you&apos;d like us to follow up on your feedback
          </p>
        </div>

        <button
          type="submit"
          disabled={status === 'submitted'}
          className="w-full bg-green-600 text-white py-2 px-4 rounded-md hover:bg-green-700 focus:outline-none focus:ring-2 focus:ring-green-500 focus:ring-offset-2 disabled:opacity-50 disabled:cursor-not-allowed transition-colors flex items-center justify-center gap-2"
        >
          {status === 'submitted' ? (
            <>
              <Check size={16} />
              Feedback Sent!
            </>
          ) : (
            <>
              <Send size={16} />
              Submit Feedback
            </>
          )}
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