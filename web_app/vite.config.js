import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { resolve } from 'path'
import fs from 'fs'

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    react(),
    {
      name: 'copy-data',
      // When running the dev server, copy the data to the public folder
      // so it can be fetched by the frontend
      configureServer(server) {
        const configPath = resolve(__dirname, 'public', '../../results/plasmid_features.json');
        const destPath = resolve(__dirname, 'public', 'plasmid_features.json');
        fs.copyFileSync(configPath, destPath);
      },
      // When building the project, copy the data to the dist folder
      // so it can be fetched by the frontend
      writeBundle() {
        const configPath = resolve(__dirname, 'public', '../../results/plasmid_features.json');
        const destPath = resolve(__dirname, 'dist', 'plasmid_features.json');
        fs.copyFileSync(configPath, destPath);
        // Write a version.env.json file with the $VERSION and $COMMIT_SHA variables
        fs.writeFileSync(resolve(__dirname, 'dist', 'version.json'), `{ "version": "${process.env.VERSION || ''}", "commit_sha": "${process.env.COMMIT_SHA || ''}" }`);
      },
    },
  ],
})
