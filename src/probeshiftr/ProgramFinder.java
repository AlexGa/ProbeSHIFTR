package probeshiftr;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class ProgramFinder {

	public String findPath(String programName) {
		
        // Check for OS and define the correct command 
		
        String os = System.getProperty("os.name").toLowerCase();
        String[] command;

        if (os.contains("win")) {
        
        	command = new String[]{"where", programName}; // Windows
        
        } else {
            
        	command = new String[]{"which", programName}; // Unix/Linux/macOS
        
        }

        try {
        	
            ProcessBuilder pb = new ProcessBuilder(command);
            Process process = pb.start();

            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String path = reader.readLine(); // Get the first line of output

            int exitCode = process.waitFor();
            
            if (exitCode == 0 && path != null && !path.isEmpty()) {
                return path.trim();
            } else {
                return null; // Program not found
            }
            
        } catch (IOException | InterruptedException e) {
            
        	e.printStackTrace();
            return null;
        
        }
    }
	
}
