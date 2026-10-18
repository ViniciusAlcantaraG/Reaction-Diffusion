#include "raylib.h"
#include <bits/stdc++.h>

using namespace std;


#define SCREEN_WIDTH 600
#define SCREEN_HEIGHT 600
#define SCALE 5

const int COL_SIZE = SCREEN_HEIGHT / SCALE;
const int ROW_SIZE = SCREEN_WIDTH / SCALE;

vector<vector<pair<float, float>>> grid(ROW_SIZE, vector<pair<float,float>>(COL_SIZE, {1.0f,0.0f}));
vector<vector<pair<float, float>>> next_grid(ROW_SIZE, vector<pair<float,float>>(COL_SIZE, {1.0f,0.0f}));

float f = 0.05;
float k = 0.06;
float D_A = 1.0;
float D_B = 0.5;

void update_reaction_diffusion();
void update_pixels(Color*);

int main() {
    InitWindow(800, 800, "Gray-Scott - raylib");
    SetTargetFPS(60);

    // 1. Create an in-memory CPU image buffer
    Image screenImage = GenImageColor(ROW_SIZE, COL_SIZE, BLACK);

    // 2. Load it onto the GPU as a texture
    Texture2D screenTexture = LoadTextureFromImage(screenImage);

    while (!WindowShouldClose()) {

        // Perform Reaction-Diffusion via Gray-Scott equations
        update_reaction_diffusion();

        // Translate the chemical concentrations into pixel color
        Color* pixels =  (Color*)screenImage.data;
        update_pixels(pixels);

        // Push CPU image bytes to the GPU texture in a single call
        UpdateTexture(screenTexture, screenImage.data);

        // Draw frame
        BeginDrawing();
            ClearBackground(BLACK);
            DrawTextureEx(screenTexture, (Vector2){0, 0}, 0.0f, 800.0f/ROW_SIZE, WHITE);
        EndDrawing();
        }

        // Step 3: Cleanup
        UnloadTexture(screenTexture);
        UnloadImage(screenImage);
        CloseWindow();

        return 0;
}

int wrap(int v, int limit){

    return ((v + limit) % limit);
}
void update_reaction_diffusion(){

    for (int r = 0; r < ROW_SIZE; r++){
        for (int c = 0; c < COL_SIZE; c++){

            next_grid[r][c] = grid[r][c];
            // Update A chemical
            // 1. Diagonal cells in the Laplacian
            next_grid[r][c].first += grid[wrap(r-1, ROW_SIZE)][wrap(c-1, COL_SIZE)].first * D_A * 0.05;
            next_grid[r][c].first += grid[wrap(r-1, ROW_SIZE)][wrap(c+1, COL_SIZE)].first * D_A * 0.05;
            next_grid[r][c].first += grid[wrap(r+1, ROW_SIZE)][wrap(c-1, COL_SIZE)].first * D_A * 0.05;
            next_grid[r][c].first += grid[wrap(r+1, ROW_SIZE)][wrap(c+1, COL_SIZE)].first * D_A * 0.05;
            
            // 2. Left, Right, Up, Down updates
            next_grid[r][c].first += grid[r][wrap(c-1, COL_SIZE)].first * D_A * 0.2;
            next_grid[r][c].first += grid[r][wrap(c+1, COL_SIZE)].first * D_A * 0.2;
            next_grid[r][c].first += grid[wrap(r-1, ROW_SIZE)][c].first * D_A * 0.2;
            next_grid[r][c].first += grid[wrap(r+1, ROW_SIZE)][c].first * D_A * 0.2; 
            
            // 3. Remove from the center cell
            next_grid[r][c].first += grid[r][c].first * D_A * -1;

            // 4. Reaction term
            next_grid[r][c].first -= grid[r][c].first * pow(grid[r][c].second, 2);

            // 5. Feed term
            next_grid[r][c].first += f*(1-grid[r][c].first);

            // Update B chemical with the same steps (Reaction is positive and kill instead of feed)
            next_grid[r][c].second += grid[wrap(r-1, ROW_SIZE)][wrap(c-1, COL_SIZE)].second * D_B * 0.05;
            next_grid[r][c].second += grid[wrap(r-1, ROW_SIZE)][wrap(c+1, COL_SIZE)].second * D_B * 0.05;
            next_grid[r][c].second += grid[wrap(r+1, ROW_SIZE)][wrap(c-1, COL_SIZE)].second * D_B * 0.05;
            next_grid[r][c].second += grid[wrap(r+1, ROW_SIZE)][wrap(c+1, COL_SIZE)].second * D_B * 0.05;
            
            next_grid[r][c].second += grid[r][wrap(c-1, COL_SIZE)].second * D_B * 0.2;
            next_grid[r][c].second += grid[r][wrap(c+1, COL_SIZE)].second * D_B * 0.2;
            next_grid[r][c].second += grid[wrap(r-1, ROW_SIZE)][c].second * D_B * 0.2;
            next_grid[r][c].second += grid[wrap(r+1, ROW_SIZE)][c].second * D_B * 0.2; 
            
            next_grid[r][c].second += grid[r][c].second * D_B * -1;

            next_grid[r][c].second += grid[r][c].first * pow(grid[r][c].second, 2);

            // Kill term
            next_grid[r][c].second -= (f+k) * grid[r][c].second;
        }
    }
    swap(grid, next_grid);
}

void update_pixels(Color* pixels){

    for (int r = 0; r < ROW_SIZE; r++){
        for (int c = 0; c < COL_SIZE; c++){
            
            float A_concentration = grid[r][c].first;
            float B_concentration = grid[r][c].second;
            
            unsigned char A_color = (unsigned char) (A_concentration * 255);
            unsigned char B_color = (unsigned char) (B_concentration * 255);

            pixels[r * COL_SIZE + c] = (Color){B_color, B_color, B_color, 255};
        }
    }
}
