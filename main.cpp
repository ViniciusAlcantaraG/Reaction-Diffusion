#include <raylib.h>
#include <bits/stdc++.h>

using namespace std;


#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 1000
#define SCALE 5

constexpr int COLS = SCREEN_WIDTH / SCALE;
constexpr int ROWS = SCREEN_HEIGHT / SCALE;
constexpr int N = COLS * ROWS;

// Chemical concentrations stored as flat array so memory access is contiguous
vector<float> A(N, 1.0f);
vector<float> B(N, 0.0f);
vector<float> A_next(N, 0.0f);
vector<float> B_next(N, 0.0f);

// CPU pixel buffer pushed to the GPU once per frame.
vector<Color> frame_buffer(N, BLACK);

// Gray-Scott feed/kill rates - settable at runtime via CLI (--f / --k),
// defaults reproduce the interactive behaviour.
float f = 0.098f;
float k = 0.057f;
const float D_A = 1.0f;
const float D_B = 0.5f;

// Time-step
constexpr float DT = 1.0f;

// B peaks around 0.3-0.4 while A sits near 1.0, so B is scaled up to be able
// to reach full yellow
constexpr float B_GAIN = 3.0f;

int steps_per_frame = 4;

void update_reaction_diffusion();
void update_pixels(Color*);
void init_grid(bool noise);

int wrap(int v, int limit){

    return ((v + limit) % limit);
}

// 9-point Laplacian: diagonals weigh 0.05, orthogonal neighbours 0.20,
// centre -1.0, so the weights sum to zero.
inline float laplacian(float centre,
                       float ul, float u, float ur,
                       float l,            float r,
                       float dl, float d, float dr){

    return 0.05f * (ul + ur + dl + dr) + 0.20f * (u + d + l + r) - centre;
}

// Reaction + feed/kill terms, shared by both chemicals (ab2 is computed once).
inline void react(float a, float b, float lap_a, float lap_b,
                  float& out_a, float& out_b){

    const float ab2 = a * b * b;
    out_a = a + DT * (D_A * lap_a - ab2 + f * (1.0f - a));
    out_b = b + DT * (D_B * lap_b + ab2 - (f + k) * b);
}

int main(int argc, char** argv) {

    // ---- CLI: optional headless export mode (no window) ----
    bool headless = false;
    bool noise_seed = false;
    int steps = 40000;
    const char* out_path = "output.png";

    for (int i = 1; i < argc; i++){
        const string arg = argv[i];
        auto next = [&](const char* what) -> const char* {
            if (i + 1 >= argc){ fprintf(stderr, "Missing value for %s\n", what); exit(1); }
            return argv[++i];
        };
        if      (arg == "--headless") headless = true;
        else if (arg == "--f")        f = stof(next("--f"));
        else if (arg == "--k")        k = stof(next("--k"));
        else if (arg == "--steps")    steps = stoi(next("--steps"));
        else if (arg == "--out")      out_path = next("--out");
        else if (arg == "--seed")     noise_seed = (string(next("--seed")) == "noise");
        else {
            fprintf(stderr, "Unknown option: %s\n", arg.c_str());
            fprintf(stderr, "Usage: %s [--headless] [--f <feed>] [--k <kill>] [--steps <n>] [--seed square|noise] [--out <file.png>]\n", argv[0]);
            return 1;
        }
    }

    if (headless){
        init_grid(noise_seed);

        for (int s = 0; s < steps; s++) update_reaction_diffusion();

        update_pixels(frame_buffer.data());

        // ExportImage works without a window: wrap the CPU pixel buffer in an
        // Image (it points at our data, no copy) and write the PNG directly.
        Image img = { frame_buffer.data(), COLS, ROWS, 1, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8 };
        ExportImage(img, out_path);
        printf("Exported %d steps (f=%.4f k=%.4f) to %s\n", steps, f, k, out_path);
        return 0;
    }

    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Gray-Scott - raylib");
    SetTargetFPS(60);

    // Create an in-memory image, upload it as a texture, then free the CPU copy
    // (the GPU owns its own data from here on).
    Image screenImage = GenImageColor(COLS, ROWS, BLACK);
    Texture2D screenTexture = LoadTextureFromImage(screenImage);
    UnloadImage(screenImage);

    // Seed the chemical fields
    init_grid(false);

    while (!WindowShouldClose()) {

        // Sub-steps per frame: one Gray-Scott step is sub-millisecond now,
        // so run several to make the pattern evolve visibly faster.
        if (IsKeyPressed(KEY_UP)   || IsKeyPressed(KEY_EQUAL)) steps_per_frame++;
        if (IsKeyPressed(KEY_DOWN) || IsKeyPressed(KEY_MINUS)) steps_per_frame = max(1, steps_per_frame - 1);

        for (int s = 0; s < steps_per_frame; s++){
            update_reaction_diffusion();
        }

        // Translate the chemical concentrations into pixel color
        update_pixels(frame_buffer.data());

        // Push CPU image bytes to the GPU texture in a single call
        UpdateTexture(screenTexture, frame_buffer.data());

        // Draw frame
        BeginDrawing();
            ClearBackground(BLACK);
            DrawTextureEx(screenTexture, (Vector2){0, 0}, 0.0f, (float)SCREEN_WIDTH / COLS, WHITE);
            DrawText(TextFormat("%d steps/frame  [UP/DOWN]", steps_per_frame), 10, 10, 20, RAYWHITE);
            DrawFPS(10, 35);
        EndDrawing();
    }

    // Cleanup
    UnloadTexture(screenTexture);
    CloseWindow();

    return 0;
}

// Seed the chemical fields: A everywhere, then a patch of B in the middle
// ("square") or scattered random patches ("noise").
void init_grid(bool noise){

    fill(A.begin(), A.end(), 1.0f);
    fill(B.begin(), B.end(), 0.0f);

    if (noise){
        // Deterministic seed so headless runs are reproducible
        mt19937 rng(1234);
        uniform_real_distribution<float> dist(0.0f, 1.0f);
        for (int r = 0; r < ROWS; r++){
            for (int c = 0; c < COLS; c++){
                if (dist(rng) < 0.02f){
                    A[r * COLS + c] = 0.0f;
                    B[r * COLS + c] = 1.0f;
                }
            }
        }
    } else {
        for (int r = (ROWS/2)-10; r < (ROWS/2)+10; r++){
            for (int c = (COLS/2)-10; c < (COLS/2)+10; c++){
                const int i = r * COLS + c;
                A[i] = 0.0f;
                B[i] = 1.0f;
            }
        }
    }
}

void update_reaction_diffusion(){

    // 1. Interior: no wrap-around needed, plain pointer arithmetic.
    //    This is the hot path - contiguous reads, no branches, vectorisable.
    for (int r = 1; r < ROWS - 1; r++){
        const int base = r * COLS;
        for (int c = 1; c < COLS - 1; c++){

            const int i = base + c;
            const float a = A[i];
            const float b = B[i];

            const float lap_a = laplacian(a,
                                          A[i - COLS - 1], A[i - COLS], A[i - COLS + 1],
                                          A[i - 1],                    A[i + 1],
                                          A[i + COLS - 1], A[i + COLS], A[i + COLS + 1]);

            const float lap_b = laplacian(b,
                                          B[i - COLS - 1], B[i - COLS], B[i - COLS + 1],
                                          B[i - 1],                    B[i + 1],
                                          B[i + COLS - 1], B[i + COLS], B[i + COLS + 1]);

            react(a, b, lap_a, lap_b, A_next[i], B_next[i]);
        }
    }

    // 2. Borders (top/bottom rows and left/right columns, including corners):
    //    the only cells that actually need the wrap-around.
    for (int r = 0; r < ROWS; r++){
        for (int c = 0; c < COLS; c++){

            if (r != 0 && r != ROWS - 1 && c != 0 && c != COLS - 1) continue;

            const int i = r * COLS + c;
            const int up_base   = wrap(r - 1, ROWS) * COLS;
            const int down_base = wrap(r + 1, ROWS) * COLS;
            const int left_col  = wrap(c - 1, COLS);
            const int right_col = wrap(c + 1, COLS);
            const int self_row  = r * COLS;

            const float a = A[i];
            const float b = B[i];

            const float lap_a = laplacian(a,
                                          A[up_base   + left_col], A[up_base   + c], A[up_base   + right_col],
                                          A[self_row + left_col],                    A[self_row + right_col],
                                          A[down_base + left_col], A[down_base + c], A[down_base + right_col]);

            const float lap_b = laplacian(b,
                                          B[up_base   + left_col], B[up_base   + c], B[up_base   + right_col],
                                          B[self_row + left_col],                    B[self_row + right_col],
                                          B[down_base + left_col], B[down_base + c], B[down_base + right_col]);

            react(a, b, lap_a, lap_b, A_next[i], B_next[i]);
        }
    }

    swap(A, A_next);
    swap(B, B_next);
}

void update_pixels(Color* pixels){

    // Chemical A -> blue (background where A is abundant),
    // Chemical B -> yellow (structures where B accumulates).
    for (int i = 0; i < N; i++){

        const float a = clamp(A[i], 0.0f, 1.0f);
        const float b = clamp(B[i] * B_GAIN, 0.0f, 1.0f);

        pixels[i] = (Color){
            (unsigned char)(b * 255.0f),
            (unsigned char)(b * 255.0f),
            (unsigned char)(a * 255.0f),
            255
        };
    }
}
