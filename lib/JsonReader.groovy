import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonParseException
import com.fasterxml.jackson.databind.JsonMappingException
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper

class JsonReader {
    static stream(String filePath, ObjectMapper mapper, Closure closure) {
        try {
            JsonFactory factory = mapper.getFactory()
            JsonParser parser = factory.createParser(new File(filePath))
            if (parser.nextToken() == JsonToken.START_ARRAY) {
                while (parser.nextToken() != JsonToken.END_ARRAY) {
                    JsonNode node = mapper.readTree(parser)
                    closure.call(node)
                }
            }
            parser.close()
        } catch (FileNotFoundException e) {
            throw new Exception("File not found: $jsonFile", e)
        } catch (JsonParseException e) {
            throw new Exception("Error parsing JSON file: $jsonFile", e)
        } catch (JsonMappingException e) {
            throw new Exception("Error mapping JSON content for file: $jsonFile", e)
        } catch (IOException e) {
            throw new Exception("IO error reading file: $jsonFile", e)
        }
    }
}
