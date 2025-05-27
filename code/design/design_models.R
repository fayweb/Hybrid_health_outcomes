# Create publication-ready table of model specifications
library(knitr)
library(kableExtra)
library(dplyr)

# Create the model specification table (updated with your actual models)
model_table <- data.frame(
    Model = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"),
    
    Name = c("Complete model", "Host factors model", "Immune signature model", 
             "Core predictors model", "Parsimonious model", "Interactive model"),
    
    PC1 = c("✓", "✗", "✓", "✓", "✓", "✓"),
    PC2 = c("✓", "✗", "✓", "✓", "✓", "✓"),
    
    `Infection Status` = c("✓", "✓", "✗", "✓", "✓", "✓"),
    `Infection Intensity` = c("✓", "✓", "✗", "✓", "✓", "✗"),
    
    `Mouse Strain` = c("✓", "✓", "✗", "✗", "✗", "✗"),
    `Immunization Status` = c("✓", "✓", "✗", "✗", "✗", "✗"),
    `Initial Body Weight` = c("✓", "✓", "✗", "✓", "✗", "✗"),
    
    `PC × Infection Interactions` = c("✗", "✗", "✗", "✗", "✗", "✓"),
    
    `R²` = c("0.631", "0.611", "0.106", "0.550", "0.547", "0.426"),
    
    Purpose = c(
        "Test all available predictors",
        "Assess non-immune factors only", 
        "Evaluate immune patterns alone",
        "Essential predictors without host factors",
        "Most parsimonious effective model",
        "Test immune-infection interactions"
    )
)

# Method 1: Basic publication table
basic_table <- model_table %>%
    kable(caption = "Model specifications for predicting maximum weight loss during Eimeria infection",
          col.names = c("Model", "Model Name", "PC1", "PC2", "Infection Status", 
                        "Infection Intensity", "Mouse Strain", "Immunization Status", 
                        "Initial Body Weight", "PC × Infection", "R²", "Purpose"),
          align = c("l", "l", rep("c", 9), "l")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE, font_size = 12) %>%
    add_header_above(c(" " = 2, "Immune Components" = 2, "Infection Variables" = 2, 
                       "Host Factors" = 3, "Interactions" = 1, "Fit" = 1, " " = 1)) %>%
    column_spec(1, bold = TRUE) %>%
    column_spec(2, italic = TRUE) %>%
    row_spec(0, bold = TRUE, background = "#f0f0f0")

print(basic_table)

# Method 2: Enhanced table with color coding (fixed footnote)
enhanced_table <- model_table %>%
    kable(caption = "Model specifications for predicting maximum weight loss during Eimeria infection",
          col.names = c("Model", "Model Name", "PC1", "PC2", "Infection Status", 
                        "Infection Intensity", "Mouse Strain", "Immunization Status", 
                        "Initial Body Weight", "PC × Infection", "R²", "Purpose"),
          align = c("l", "l", rep("c", 9), "l")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = FALSE, font_size = 11) %>%
    add_header_above(c(" " = 2, "Immune Components" = 2, "Infection Variables" = 2, 
                       "Host Factors" = 3, "Interactions" = 1, "Fit" = 1, " " = 1)) %>%
    column_spec(1, bold = TRUE, color = "black", background = "#e8f4f8") %>%
    column_spec(2, italic = TRUE, width = "2.5cm") %>%
    column_spec(3:4, background = "#f0f8e8") %>%  # Light green for immune components
    column_spec(5:6, background = "#fff8e8") %>%  # Light yellow for infection variables  
    column_spec(7:9, background = "#f8f0f8") %>%  # Light purple for host factors
    column_spec(10, background = "#e8f0f8") %>%   # Light blue for interactions
    column_spec(11, background = "#f0f0f0") %>%   # Light gray for R²
    column_spec(12, width = "3cm") %>%
    row_spec(0, bold = TRUE, background = "#d0d0d0", color = "black") %>%
    footnote(general = "✓ indicates variable included; ✗ indicates variable excluded")

print(enhanced_table)

# Method 3: Compact version for journals
compact_table <- model_table %>%
    select(-Purpose) %>%  # Remove purpose for space
    kable(caption = "Model specifications for maximum weight loss prediction",
          col.names = c("Model", "Name", "PC1", "PC2", "Infection Status", 
                        "Infection Intensity", "Mouse Strain", "Immunization", 
                        "Initial Weight", "Interactions", "R²"),
          align = c("l", "l", rep("c", 9))) %>%
    kable_styling(bootstrap_options = c("striped", "condensed"),
                  full_width = FALSE, font_size = 10) %>%
    add_header_above(c(" " = 2, "Immune" = 2, "Infection" = 2, "Host Factors" = 3, " " = 2)) %>%
    column_spec(1, bold = TRUE) %>%
    column_spec(2, italic = TRUE, width = "2cm") %>%
    footnote(general = "✓ = included; ✗ = excluded")

print(compact_table)

# Save tables
save_kable(enhanced_table, "output/tables/model_specifications_enhanced.html")
save_kable(compact_table, "output/tables/model_specifications_compact.html")

# Print summary
cat("Tables created successfully!\n")
cat("Key findings from your models:\n")
cat("- Complete model (Model 1): R² = 0.631 (best overall fit)\n")
cat("- Interactive model (Model 6): R² = 0.426 (significant interactions found)\n")
cat("- Immune signature alone (Model 3): R² = 0.106 (modest but significant)\n")
cat("- Host factors model (Model 2): R² = 0.611 (nearly as good without immune data)\n")

# Create model performance summary
model_performance <- data.frame(
    Model = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"),
    Name = c("Complete", "Host factors", "Immune signature", "Core predictors", "Parsimonious", "Interactive"),
    R2 = c(0.631, 0.611, 0.106, 0.550, 0.547, 0.426),
    Variables = c(22, 20, 2, 6, 4, 8),
    AIC_approx = c(NA, NA, NA, NA, NA, NA)  # You can add these if you calculated them
)

# Save performance summary
write.csv(model_performance, "output/tables/model_performance_summary.csv", row.names = FALSE)

