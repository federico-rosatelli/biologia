import {fileURLToPath, URL} from 'node:url'

import {defineConfig} from 'vite'
import vue from '@vitejs/plugin-vue'

// https://vitejs.dev/config/
export default defineConfig(({command, mode, ssrBuild}) => {
	const ret = {
		plugins: [vue()],
		resolve: {
			alias: {
				'@': fileURLToPath(new URL('./src', import.meta.url))
			}
		},
		devServer: {
			proxy: 'http://localhost:8080/',
			port:8080,
		},
		publicPath: JSON.stringify("http://localhost:8080")
	};
	if (command === 'serve') {
		ret.define = {
			"__API_URL__": JSON.stringify("http://localhost:8080"),
		};
	} else {
		ret.define = {
			"__API_URL__": JSON.stringify("http://localhost:8080"),
		};
	}
	return ret;
})