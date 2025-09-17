#!/bin/bash

echo "🚀 Setting up GitHub repository for automatic Pages deployment"
echo ""

# Check if we're in a git repository
if [ ! -d ".git" ]; then
    echo "📁 Initializing Git repository..."
    git init
    echo "✅ Git repository initialized"
else
    echo "✅ Git repository found"
fi

# Add all files
echo "📝 Adding files to Git..."
git add .
git commit -m "Initial commit: Bayesian SIR WebAssembly application with GitHub Pages deployment"

echo ""
echo "🔧 Next steps to deploy to GitHub Pages:"
echo ""
echo "1. Create a new repository on GitHub:"
echo "   https://github.com/new"
echo ""
echo "2. Replace YOUR_USERNAME in README.md with your GitHub username"
echo ""
echo "3. Add your GitHub repository as remote:"
echo "   git remote add origin https://github.com/YOUR_USERNAME/sir_bayes.git"
echo ""
echo "4. Push to GitHub:"
echo "   git push -u origin main"
echo ""
echo "5. Enable GitHub Pages in repository settings:"
echo "   - Go to Settings → Pages"
echo "   - Source: GitHub Actions"
echo "   - The workflow will automatically build and deploy"
echo ""
echo "6. Your app will be live at:"
echo "   https://YOUR_USERNAME.github.io/sir_bayes/"
echo ""
echo "🎉 The GitHub Actions workflow will automatically:"
echo "   ✅ Build the C++ WebAssembly module"
echo "   ✅ Deploy to GitHub Pages"
echo "   ✅ Update on every push to main branch"
