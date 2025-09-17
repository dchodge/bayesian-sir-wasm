// src/server.cpp - Simple C++ HTTP Server for WebR Application
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <thread>
#include <filesystem>
#include <cstring>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <signal.h>

class HttpServer {
private:
    int server_fd;
    int port;
    std::string web_root;
    std::map<std::string, std::string> mime_types;
    bool running;

    void initializeMimeTypes() {
        mime_types[".html"] = "text/html";
        mime_types[".js"] = "application/javascript";
        mime_types[".css"] = "text/css";
        mime_types[".wasm"] = "application/wasm";
        mime_types[".json"] = "application/json";
        mime_types[".png"] = "image/png";
        mime_types[".jpg"] = "image/jpeg";
        mime_types[".gif"] = "image/gif";
        mime_types[".ico"] = "image/x-icon";
    }

    std::string getMimeType(const std::string& filepath) {
        size_t dot = filepath.find_last_of('.');
        if (dot != std::string::npos) {
            std::string ext = filepath.substr(dot);
            auto it = mime_types.find(ext);
            if (it != mime_types.end()) {
                return it->second;
            }
        }
        return "text/plain";
    }

    std::string readFile(const std::string& filepath) {
        std::ifstream file(filepath, std::ios::binary);
        if (!file.is_open()) {
            return "";
        }
        
        std::ostringstream contents;
        contents << file.rdbuf();
        return contents.str();
    }

    std::string createResponse(int status_code, const std::string& status_text, 
                              const std::string& content_type, const std::string& body) {
        std::ostringstream response;
        response << "HTTP/1.1 " << status_code << " " << status_text << "\r\n";
        response << "Content-Type: " << content_type << "\r\n";
        response << "Content-Length: " << body.length() << "\r\n";
        
        // CORS headers for WebAssembly
        response << "Cross-Origin-Embedder-Policy: require-corp\r\n";
        response << "Cross-Origin-Opener-Policy: same-origin\r\n";
        response << "Cross-Origin-Resource-Policy: cross-origin\r\n";
        response << "Access-Control-Allow-Origin: *\r\n";
        response << "Access-Control-Allow-Methods: GET, POST, OPTIONS\r\n";
        response << "Access-Control-Allow-Headers: Content-Type\r\n";
        
        // Content Security Policy for WebR (allows inline scripts and eval for development)
        response << "Content-Security-Policy: default-src 'self'; script-src 'self' 'unsafe-eval' 'unsafe-inline'; worker-src 'self' blob:; style-src 'self' 'unsafe-inline'; connect-src 'self'; img-src 'self' data:; font-src 'self' data: https: moz-extension: chrome-extension:;\r\n";
        
        response << "\r\n";
        response << body;
        
        return response.str();
    }

    void handleClient(int client_socket) {
        char buffer[4096] = {0};
        ssize_t bytes_read = read(client_socket, buffer, sizeof(buffer) - 1);
        
        if (bytes_read <= 0) {
            close(client_socket);
            return;
        }

        std::string request(buffer, bytes_read);
        std::istringstream request_stream(request);
        std::string method, path, version;
        request_stream >> method >> path >> version;

        // Handle default path
        if (path == "/") {
            path = "/index.html";
        }

        // Security: prevent directory traversal
        if (path.find("..") != std::string::npos) {
            std::string response = createResponse(403, "Forbidden", "text/plain", "Access denied");
            send(client_socket, response.c_str(), response.length(), 0);
            close(client_socket);
            return;
        }

        std::string full_path = web_root + path;
        
        // Check if file exists
        if (!std::filesystem::exists(full_path)) {
            std::string response = createResponse(404, "Not Found", "text/html", 
                "<html><body><h1>404 Not Found</h1><p>The requested file was not found.</p></body></html>");
            send(client_socket, response.c_str(), response.length(), 0);
            close(client_socket);
            return;
        }

        // Read and serve file
        std::string content = readFile(full_path);
        if (content.empty()) {
            std::string response = createResponse(500, "Internal Server Error", "text/plain", "Error reading file");
            send(client_socket, response.c_str(), response.length(), 0);
            close(client_socket);
            return;
        }

        std::string mime_type = getMimeType(full_path);
        std::string response = createResponse(200, "OK", mime_type, content);
        
        send(client_socket, response.c_str(), response.length(), 0);
        close(client_socket);
    }

public:
    HttpServer(int port, const std::string& web_root) : port(port), web_root(web_root), running(false) {
        initializeMimeTypes();
    }

    bool start() {
        // Create socket
        server_fd = socket(AF_INET, SOCK_STREAM, 0);
        if (server_fd == -1) {
            std::cerr << "Failed to create socket\n";
            return false;
        }

        // Set socket options
        int opt = 1;
        if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt))) {
            std::cerr << "Failed to set socket options\n";
            close(server_fd);
            return false;
        }

        // Bind socket
        struct sockaddr_in address;
        address.sin_family = AF_INET;
        address.sin_addr.s_addr = INADDR_ANY;
        address.sin_port = htons(port);

        if (bind(server_fd, (struct sockaddr*)&address, sizeof(address)) < 0) {
            std::cerr << "Failed to bind socket to port " << port << "\n";
            close(server_fd);
            return false;
        }

        // Listen for connections
        if (listen(server_fd, 10) < 0) {
            std::cerr << "Failed to listen on socket\n";
            close(server_fd);
            return false;
        }

        running = true;
        std::cout << "🌊 WebR Oscillating Sine Wave Server\n";
        std::cout << "✅ Server started on http://localhost:" << port << "\n";
        std::cout << "📁 Serving from: " << std::filesystem::absolute(web_root) << "\n";
        std::cout << "🔄 WebR + C++ + WebAssembly\n";
        std::cout << "Press Ctrl+C to stop\n\n";

        return true;
    }

    void run() {
        while (running) {
            struct sockaddr_in client_addr;
            socklen_t client_len = sizeof(client_addr);
            
            int client_socket = accept(server_fd, (struct sockaddr*)&client_addr, &client_len);
            if (client_socket < 0) {
                if (running) {
                    std::cerr << "Failed to accept client connection\n";
                }
                continue;
            }

            // Handle client in a separate thread for better performance
            std::thread client_thread(&HttpServer::handleClient, this, client_socket);
            client_thread.detach();
        }
    }

    void stop() {
        running = false;
        if (server_fd != -1) {
            close(server_fd);
        }
    }
};

// Global server instance for signal handling
HttpServer* global_server = nullptr;

void signal_handler(int signal) {
    std::cout << "\n🛑 Server stopping...\n";
    if (global_server) {
        global_server->stop();
    }
    exit(0);
}

int main(int argc, char* argv[]) {
    int port = 8080;
    std::string web_root = "web";

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--port" && i + 1 < argc) {
            port = std::stoi(argv[i + 1]);
            i++;
        } else if (std::string(argv[i]) == "--dir" && i + 1 < argc) {
            web_root = argv[i + 1];
            i++;
        } else if (std::string(argv[i]) == "--help") {
            std::cout << "Usage: " << argv[0] << " [options]\n";
            std::cout << "Options:\n";
            std::cout << "  --port <port>    Set server port (default: 8080)\n";
            std::cout << "  --dir <path>     Set web root directory (default: web)\n";
            std::cout << "  --help           Show this help message\n";
            return 0;
        }
    }

    // Verify web directory exists
    if (!std::filesystem::exists(web_root)) {
        std::cerr << "❌ Error: Web directory '" << web_root << "' not found!\n";
        std::cerr << "Please make sure the web directory exists with your HTML files.\n";
        return 1;
    }

    // Set up signal handling
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    // Create and start server
    HttpServer server(port, web_root);
    global_server = &server;

    if (!server.start()) {
        std::cerr << "❌ Failed to start server\n";
        return 1;
    }

    // Open browser
    std::string open_cmd;
#ifdef __APPLE__
    open_cmd = "open http://localhost:" + std::to_string(port);
#elif __linux__
    open_cmd = "xdg-open http://localhost:" + std::to_string(port);
#else
    open_cmd = "start http://localhost:" + std::to_string(port);
#endif
    
    std::cout << "🚀 Opening browser...\n";
    system(open_cmd.c_str());

    server.run();
    return 0;
}
