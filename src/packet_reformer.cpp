#include "packet_reformer.h"
#include "crcinterface.h"
#include <bit>

PacketReformer::PacketReformer(
    std::vector<bool> * in_data_buffer,
    std::vector<SimpPacket> * control_packet_buffer,
    std::vector<SimpPacket> * data_packet_buffer) : 
    inDataBuffer(in_data_buffer),
    controlPacketBuffer(control_packet_buffer),
    dataPacketBuffer(data_packet_buffer) {

    // Make crc generator with polynomial for ISO 3309 in reversed format
    crcGenny = crcutil_interface::CRC::Create(0xEDB88320, 0, 32, true, 0, 0, 0, true, NULL);

}

std::vector<bool>::iterator PacketReformer::getInputPosition() {

    return inputPos;

}

uint8_t PacketReformer::formPacket(bool are_we_source, uint8_t sender_address, uint8_t receiver_address) {

    // Temp data storage
    uint8_t tempFlag;
    uint8_t tempLen;
    bool noFlag;
    SimpPacket tempPacket;
    uint8_t tempByte;
    uint8_t tempRecAddress;
    uint8_t tempSendAddress;
    std::vector<bool>::iterator endPos;

    // CRC
    crcutil_interface::UINT64 tempCRC;
    uint32_t erc;

    // Clear temporary packet and other data
    (tempPacket.data).clear();
    tempPacket.controlCode = 0;
    tempPacket.dataLength = 0;
    tempPacket.sequenceNumber = 0;
    tempPacket.erc = false;

    tempFlag = 0;
    tempLen = 0;
    noFlag = true;
    tempByte = 0;
    tempRecAddress = 0;
    tempSendAddress = 0;

    // Setup
    noFlag = true;
    inputPos = inDataBuffer->begin();

    // Find flag
    while(noFlag) {

        // Check Buffer has enough data for flag (8 bits)
        if(std::distance(inputPos, inDataBuffer->end()) <= 8) {

            inputPos = inDataBuffer->begin();
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

    // Pointer currently one bit into flag

    // Check for second flag

    // Check Buffer has enough data for length field (8 bits for control + 7 bits for rest of flag)
    if(std::distance(inputPos + 7, inDataBuffer->end()) <= 8) {

        inputPos = inDataBuffer->begin();
        return 0;
        
    }

    // Move pointer to one bit past end of flag field
    inputPos += 7;

    // Read potential len field: inputPos one bit after flag, need to check 4 forward for length
    tempLen = 0;
    for(int i = 0; i < 4; i++) {

        tempLen += (*((inputPos+4) + i) << (3 - i));

    }

    // Second flag should be 4 control bits + 4 length bits + (8 + 8) address bits + 8 sequence bits + 8*length data bits + 32 CRC bits ahead of inputPos

    // Check that enough data is present in buffer
    if(std::distance(inputPos, inDataBuffer->end()) < 4 + 4 + 16 + 8 + 8*tempLen + 32 + 8) {

        inputPos = inDataBuffer->begin();
        return 0;
        
    }

    // We now know enough data in buffer, inputPos one bit after flag

    // Read potential flag field
    tempFlag = 0;
    for(int i = 0; i < 8; i++) {

        tempFlag += (*((inputPos + 4 + 4 + 8 + 8 + 8 + 8*tempLen + 32) + i) << (7 - i));

    }

    if(tempFlag != 0b01111110) {

        // At a minimum, our first flag is wrong, so erase to there from buffer
        inDataBuffer->erase(inDataBuffer->begin(), inputPos - 1);
        return 1;

    }

    // Now we have verified that both flags are good, so it is time to make the SimpPacket
    // We know that the packet is the correct length, no more need to check size

    // Let us store position of end of packet for later - we know all packets are 'valid' now, so we should erase to end of packet +1 from now on
    // input pos still one bit after flag
    endPos = inputPos + 4 + 4 + 8 + 8 + 8 + 8*tempLen + 32 + 8;

    // Read Control
    for(int i = 0; i < 4; i++) {

        tempPacket.controlCode += ((*inputPos) << (3 - i));
        inputPos++;

    }

    // Read Data Length - except we alread did, it is tempLen
    tempPacket.dataLength = tempLen;
    inputPos += 4;

    // Read Addresses

    // Reciever
    for(int i = 0; i < 8; i++) {

       tempRecAddress += ((*inputPos) << (8 - i));
       inputPos++;

    }

    // Check if our address, return appropriate value if not
    if(tempRecAddress != receiver_address) {

        inputPos = endPos;
        inDataBuffer->erase(inDataBuffer->begin(), endPos);
        return 2;

    }

    // Sender

    for(int i = 0; i < 8; i++) {

        tempSendAddress += ((*inputPos) << (8 - i));
        inputPos++;

    }

    // Check if sender we expect, return appropriate value if not
    if(tempSendAddress != sender_address) {

        inputPos = endPos;
        inDataBuffer->erase(inDataBuffer->begin(), endPos);
        return 3;
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

    // If data length isn't zero, check checksum
    if(tempPacket.dataLength != 0) {

        // Read Checksum if data length != 0
        erc = 0;
        for(int j = 0; j < 32; j++) {

            erc += ((*inputPos) << (32 - j));
            inputPos++;

        }

        // Calculate Checksum
        crcGenny->Compute(tempPacket.data.data(), tempPacket.dataLength, &tempCRC, NULL);
        
        // See if calculated and sent checksum match
        tempPacket.erc = (erc == (uint32_t)tempCRC);

    } else { // Otherwise act like good checksum, move pointer past checksum field

        tempPacket.erc = true;
        inputPos += 32;

    }

    // Move pointer past flag
    inputPos += 8;

    // Put to appropriate buffer
    sortPacket(tempPacket, are_we_source);
    
    // Clear the input data buffer to appropriate point
    inputPos = endPos;    
    inDataBuffer->erase(inDataBuffer->begin(), endPos);

    // Return appropriate value for good checksum else return value for bad checksum
    if(!tempPacket.erc) {

        return 255;

    } else {

        return 254;
        
    }

}

void PacketReformer::sortPacket(SimpPacket to_sort, bool are_we_source) {

    // If we are sink and code 0000, then data packet, else control
    if((to_sort.controlCode == 0b0000) && !are_we_source) {

        dataPacketBuffer->push_back(to_sort);

    } else {

        controlPacketBuffer->push_back(to_sort);

    }

}