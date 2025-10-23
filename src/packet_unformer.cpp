#include "packet_unformer.h"
#include "crcinterface.h"
#include <bit>

PacketUnformer::PacketUnformer(
    std::vector<bool> * in_data_buffer,
    std::vector<uint8_t> *  out_data_buffer) : 
    inDataBuffer(in_data_buffer),
    outDataBuffer(out_data_buffer) {

    // Make crc generator with polynomial for ISO 3309 in reversed format
    crcGenny = crcutil_interface::CRC::Create(0xEDB88320, 0, 32, true, 0, 0, 0, true, NULL);

}

uint8_t PacketUnformer::formPacket() {

    // Temp data storage
    uint8_t tempFlag;
    uint8_t tempLen;
    bool noFlag;
    std::vector<bool>::iterator inputPos;
    SimpPacket tempPacket;
    uint8_t tempByte;

    // CRC
    crcutil_interface::UINT64 tempCRC;
    uint32_t erc;

    // Clear temporary packet
    (tempPacket.data).clear();
    tempPacket.controlCode = 0;
    tempPacket.dataLength = 0;
    tempPacket.recieverAddress = 0;
    tempPacket.senderAddress = 0;
    tempPacket.sequenceNumber = 0;
    tempPacket.erc = false;

    // Setup
    noFlag = true;
    inputPos = inDataBuffer->begin();

    // Find flag
    while(noFlag) {

        // Check Buffer has enough data for flag (8 bits)
        if(std::distance(inputPos, inDataBuffer->end()) <= 8) {

            return 0;
        }

        // First read potential flag field
        tempFlag = 0;
        for(int i = 0; i < 8; i++) {

            tempFlag += (*(inputPos + i) << (7 - i));

        }

        // Check if candidate flag good
        noFlag = (tempFlag != 0b01111110);

        // Move up pointer
        inputPos++;

    }

    // Move pointer to one bit past end of flag field
    inputPos += 7;

    // Check for second flag

    // Check Buffer has enough data for length field (8 bits)
    if(std::distance(inputPos, inDataBuffer->end()) <= 8) {

        return 0;
    }

    // Read potential len field: inputPos one bit after flag, need to check 4 forward for length
    tempLen = 0;
    for(int i = 0; i < 4; i++) {

        tempLen += (*((inputPos+4) + i) << (3 - i));

    }

    // Second flag should be 4 control bits + 4 length bits + (8 + 8) address bits + 8 sequence bits + 8*length data bits + 32 CRC bits ahead of inputPos

    // Check that enough data is present in buffer
    if(std::distance(inputPos, inDataBuffer->end()) < 4 + 4 + 16 + 8 + 8*tempLen + 32 + 8) {

        return 0;
        
    }
    

    // Read potential flag field
    tempFlag = 0;
    for(int i = 0; i < 8; i++) {

        tempFlag += (*((inputPos + 4 + 4 + 8 + 8 + 8 + 8*tempLen + 32) + i) << (7 - i));

    }

    if(tempFlag != 0b01111110) {

        return 1;

    }

    // Now we have verified that both flags are good, so it is time to make the SimpPacket
    // We know that the packet is the correct length, no more need to check size

    // Read Control
    for(int i = 0; i < 4; i++) {

        tempPacket.controlCode += ((*inputPos) << (3 - i));
        inputPos++;

    }

    // Read Data Length
    for(int i = 0; i < 4; i++) {

        tempPacket.dataLength += ((*inputPos) << (3 - i));
        inputPos++;

    }

    // Read Addresses

    // Reciever
    for(int i = 0; i < 8; i++) {

        tempPacket.recieverAddress += ((*inputPos) << (8 - i));
        inputPos++;

    }

    // Sender
    for(int i = 0; i < 8; i++) {

        tempPacket.senderAddress += ((*inputPos) << (8 - i));
        inputPos++;

    }

    // Read Sequence Number
    for(int i = 0; i < 8; i++) {

        tempPacket.sequenceNumber += ((*inputPos) << (8 - i));
        inputPos++;

    }

    // Read Data

    // Outer loop is bytes, inner loop is bits
    for(int i = 0; i < tempPacket.dataLength; i++) {

        // Read next byte
        tempByte = 0;
        for(int j = 0; j < 8; j++) {

            tempByte += ((*inputPos) << (8 - j));
            inputPos++;

        }

        (tempPacket.data).push_back(tempByte);

    }

    // Read Checksum
    erc = 0;
    for(int j = 0; j < 32; j++) {

        erc += ((*inputPos) << (32 - j));
        inputPos++;

    }

    // Calculate Checksum
    crcGenny->Compute(tempPacket.data.data(), tempPacket.dataLength, &tempCRC, NULL);
    
    // See if calculated and sent checksum match
    tempPacket.erc = (erc == (uint32_t)tempCRC);

    // Return appropriate value if bad checksum, else return appropriate value for good checksum
    if(!tempPacket.erc) {

        return 2;

    } else {

        return 255;
        
    }

}