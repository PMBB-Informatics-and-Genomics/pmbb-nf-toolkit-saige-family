# pmbb-nf-toolkit-saige-family


#### The report HTML File which displays the results for each EXWAS will run as follow:
1. General Steps to follow if you are developing report on your data:
   
   a. Create a folder for keeping all the files along the html page.
   
   b. Store the Plots generated inside the Plots folder in the following format and keep the plots folder inside the folder created in step a.
   c. Generate a Manifest file plots_manifest.csv file in following format and keep it inside the folder created 
   d. create the manifest file for all analysis types in a single file (saige_exwas_suggestive_regions.csv) and keep in the folder created in step a. 
   
3. Step 2 is to run the following command to display the results. 
         
```
  1. First run the command in terminal in the directory of exwas_html_report
     python3 -m http.server 8003
  2. The file and results will be available on
     http://localhost:8003/Nextflow-static-results.html
```

#### EXWAS Report Outlook

<img width="1046" alt="Screenshot 2024-09-23 at 1 55 35 PM" src="https://github.com/user-attachments/assets/ffc3a3b4-8100-4f44-bd5c-f4b6180b11ed">

<img width="1017" alt="Screenshot 2024-09-23 at 1 55 55 PM" src="https://github.com/user-attachments/assets/3145c543-9c52-49f8-a879-1f7ec00560c5">

<img width="1037" alt="Screenshot 2024-09-23 at 1 56 05 PM" src="https://github.com/user-attachments/assets/fe9a0167-65a0-4181-9c8b-c079a5c319f6">
